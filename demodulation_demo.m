%% Passes the audio signal through a filter bank and then demodulates the subbands using GPPAD to calculate the amplitude envelopes

addpath(genpath('sinemodel'))
addpath(genpath('filterbank'))
addpath(genpath('gppad'))

fileName = 'stim7_applause'; % name of wav file
suffix = ''; % suffix for .mat file to be saved (in case of multiple versions)
numFilt = 16; % number of subbands (min. 14, but don't set too high)
len = [250 250]; % range of demodulation lengthscales (different frequencies tend to need different lengthscales)
loadPrev = 0; % load previously saved demodulation results? 
hil4lo = 0; % use the Hilbert envelope for this many subbands (starting from lowest)
filtType = 'ERB'; % 'linear' or 'ERB'

if loadPrev == 1
    matNameDemod = strcat('audio_lfm_data/demod_',fileName,suffix,'.mat'); % name of demodulation .mat file
    demod_struct = load(matNameDemod);
    % Plot
    figure(1);clf
    subplot(211)
    plot(demod_struct.A)
    title('GP env.')
    subplot(212)
    plot(demod_struct.A_)
    title('Hilbert env.')
else
    demod_struct = demod_filterbank(fileName,suffix,numFilt,loadPrev,len,hil4lo,filtType);
    save(strcat('audio_lfm_data/demod_',fileName,suffix,'.mat'),'-struct','demod_struct');
end

%% Check NMF resynthesis
err = [];
err_diff = [];
recon_sig = [];
for numb=1:8
    [w,h] = nnmf(demod_struct.A,numb);
    recon = w*h;
    err = [err; sum(sum((demod_struct.A - recon).^2)) / size(demod_struct.A,1)]; %#ok<AGROW>
    if numb>1
        err_diff = [err_diff; err(end)-err(end-1)]; %#ok<AGROW>
    else
        err_diff = 0;
    end
    recon_sig = [recon_sig; sum(recon'.*demod_struct.carrier_synth')]; %#ok<AGROW>
end
%%
sound(recon_sig(3,:),demod_struct.fs)

%% Synthesis
pind = 1:demod_struct.D;
Y = demod_struct.Y;
A = demod_struct.A;
carrier_synth = demod_struct.carrier_synth;
fs = demod_struct.fs;
Ysynth = sum(Y(:,pind),2);
Asynth = sum(A(:,pind)'.*carrier_synth(:,pind)');
sound(Ysynth,fs)
pause(size(Y,1)/fs + 0.5)
sound(Asynth,fs)