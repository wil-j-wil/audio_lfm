function demod_struct = demod_filterbank(fileName,suffix,D,loadDemod,lenRange,hil4lo,filtType)

% performs the filtering and demodulation

% Will Wilkinson - Jan 2018

if nargin < 3 || isempty(D)
    D = 16;
end
if nargin < 4 || isempty(loadDemod)
    loadDemod = 0;
end
if nargin < 5 || isempty(lenRange)
    lenRange = [100 200];
end
if nargin < 6 || isempty(hil4lo)
    hil4lo = 0;
end
if nargin < 7 || isempty(filtType)
    filtType = 'ERB';
end


% Use GPPAD to do sub-band demodulation.
% Take a sound, filter it through a Gammatone Filter bank,
% and then use GPMAP.m to demodulate the filter bank.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Loading the sentence sound

% if loadDemod is 1 then just load a previously run file, else run GPPAD
if loadDemod
    current_dir = pwd;
    if strcmp(current_dir(end-10:end),'experiments')
        demod_struct = load(strcat('../audio_lfm_data/demod_',fileName,suffix,'.mat'));
    else
        demod_struct = load(strcat('audio_lfm_data/demod_',fileName,suffix,'.mat'));
    end
    A = demod_struct.A;
    C = demod_struct.C;
    D = demod_struct.D;
    fs = demod_struct.fs;
    A_ = demod_struct.A_;
    Fc = demod_struct.Fc;
    y = demod_struct.y;
    
else
    
    current_dir = pwd;
    if strcmp(current_dir(end-10:end),'experiments')
        [y,fs] = audioread(strcat('../audio/',fileName,'.wav')); % reads in the file
    else
        [y,fs] = audioread(strcat('audio/',fileName,'.wav')); % reads in the file
    end
    %y = y(1:20000);
    
    if fs > 24000
        DS = 2; % Down sample the original if desired
    else 
        DS = 1;
    end
    y = y(1:DS:end); 
    fs = fs/DS;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Filtering the sound

    % Parameters of the filter bank
    %FLim = [10,10000]; % upper and lower frequency limits on the filters
    %D = 16; % Number of channels in the filter bank
    
    % Construct filterbank
    %[Y,Fc,g,a,c] = ltfat_filterbank(y,fs,D+1);drawnow(); % We actually get D+2 subbands but we will discard the outer 2
    %[audio_filts, audio_cutoffs_Hz] = make_erb_cos_filters(synth_dur_smp, P.audio_sr, P.N_audio_channels, P.low_audio_f, P.hi_audio_f);
    if strcmp(filtType,'ERB')
        [audio_filts, Fc] = make_erb_cos_filters(length(y), fs, D, 20, fs/2);
    elseif strcmp(filtType,'linear')
        [audio_filts, Fc] = make_lin_cos_filters(length(y), fs, D, 20, fs/2);
    else
        error('filter type not recognised')
    end
    
    Y = generate_subbands(y, audio_filts); % audio subbands
    
    if max(max(isnan(Y))) > 0
        disp('the filtered subbands contains NaNs')
    end
    Y = Y(:,2:D+1);
    Fc = Fc(2:D+1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pick the time-scale of the demodulation
    
    %lenRange = [1000 2000];
    %lenRange = [500 1000];
    %lenRange = [50 200];
    %FcRange = [Fc(1) Fc(end)];
    %len = lenRange(end)-( (lenRange(end)-lenRange(1)) * Fc/Fc(end) );
    if length(lenRange(end):lenRange(1)) == 1
        len = lenRange(1)*ones(1,D);
    else
        if lenRange(1) < lenRange(end)
            lenRange = fliplr(lenRange);
        end
        len = round(lenRange(end):-(lenRange(end)-lenRange(1))/(D-1):lenRange(1));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Run the demodulation algorithm
    
    [A,C] = GPPAD(Y,len);
    
    A_ = abs(hilbert(Y)); % subband Hilbert envelopes
    %subband_phases = hilbert(Y)./A_;
    C_ = Y./A_;
    
    if hil4lo > 0
        A(:,1:hil4lo) = A_(:,1:hil4lo);
        C(:,1:hil4lo) = C_(:,1:hil4lo);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save results
    demod_struct = struct;
    demod_struct.fileName = fileName;
    demod_struct.A = A;
    demod_struct.C = C;
    demod_struct.D = D;
    demod_struct.audio_filts = audio_filts;
    demod_struct.Fc = Fc;
    demod_struct.DS = DS;
    demod_struct.len = len;
    demod_struct.y = y;
    demod_struct.Y = Y;
    demod_struct.fs = fs;
    demod_struct.A_ = A_;
    demod_struct.C_ = C_;
    demod_struct.hil4lo = hil4lo;

end

figure(1);clf
subplot(211)
plot(A)
title('Gaussian Process Envelopes (GPPAD)')
subplot(212)
plot(A_)
title('Hilbert Envelopes')
drawnow();

%% Analyse Subband Carriers
carriers = cell(D,1);
AR_order = 120;
envLen = max(350-floor(Fc),50); % length of smoothing/normalisation envelope
carrier_synth = zeros(size(C));
if length(y) < 5000
    fft_size = 512;
    win_size = 512;
    hop_size = 256;
elseif length(y) < 10000
    fft_size = 1024;
    win_size = 1024;
    hop_size = 512;
else
    fft_size = 2048;
    win_size = 2048;
    hop_size = 1024;
end
tic
for i=1:D
    fprintf('synthesising carrier signal %d/%d\n',i,D)
    [sinu,res,F,M] = sin_res(C(:,i),fs,fft_size,win_size,hop_size);
    %[~,mind] = max(max(M,[],2));% only use the most prominent sinusoid
    %F = F(mind,:);
    %M = M(mind,:);
    %sinf = median(F(find(isnan(F)==0))); %#ok<FNDSB>
    %sinm = median(M(find(isnan(M)==0))); %#ok<FNDSB>
    %sin_synth = sinm*sin(2*pi*sinf*((1:length(res))+rand(1))/fs)'; % random initial phase
    
    %F = F(1:3,:);
    %M = M(1:3,:);
    sinf = zeros(size(M,1),1);
    sinm = zeros(size(M,1),1);
    sin_synth = zeros(size(M,1),length(res));
    for j=1:size(M,1)
        F_ = F(j,:);
        M_ = M(j,:);
        sinf(j) = median(F_(find(isnan(F_)==0))); %#ok<FNDSB>
        sinm(j) = median(M_(find(isnan(M_)==0))); %#ok<FNDSB>
        sin_synth(j,:) = sinm(j)*sin(2*pi*sinf(j)*((1:length(res))+rand(1))/fs); % random initial phase
    end
    [~,mind]=sort(sinm,'descend');
    if isempty(mind)
        mind_ = mind;
    else
        mind_ = mind(1);
    end
    for k=2:length(mind)
        fdiff = abs(sinf(mind(k)) - sinf(mind(1:k-1)));
        if fdiff > 40 % must be 40Hz apart
            mind_ = [mind_;mind(k)]; %#ok<AGROW>
        end
    end
    
    if length(mind_) > 1
        sin_synth = sum(sin_synth(mind_(1:min(5,length(mind_))),:))';
    else
        sin_synth = sin_synth(mind_,:)';
    end
    
    % aliasing sometimes occurs in the residual - filter it out
    %g_ = cell(1,1);
    %g_{1} = g{i+1}; % it's +1 because we discard the first and last filters
    %a_ = a(i+1,:);
    %c_=filterbank(res,{'realdual',g},a);
    %c__ = cell(1,1);
    %c__{1} = c_{i+1};
    %r_=2*real(ifilterbank(c__,g_,a_));
    %res=r_(1:length(res))*2;
    
    coeff = lpc(res,AR_order);
    ar_model = arima('Constant',0,'AR',-coeff(2:end),'Variance',.00001);
    res_synth = simulate(ar_model,size(res,1));
    
    %c_=filterbank(res_synth,{'realdual',g},a);
    %c__ = cell(1,1);
    %c__{1} = c_{i+1};
    %r_=2*real(ifilterbank(c__,g_,a_));
    %res_synth=r_(1:length(res))*2;
    
    res_synth = median(envelope(res,envLen(i),'rms')) * res_synth ./ median(envelope(res_synth,envLen(i),'rms'));
    
    if isempty(sin_synth)
        sig_synth = res_synth;
    else
        sig_synth = sin_synth + res_synth;
    end
    
    carriers{i} = struct;
    carriers{i}.sig = C(:,i);
    carriers{i}.sin = sinu(1:length(y));
    carriers{i}.res = res;
    carriers{i}.sinf = sinf;
    carriers{i}.sinm = sinm;
    carriers{i}.sin_synth = sin_synth;
    carriers{i}.res_synth = res_synth;
    carriers{i}.sig_synth = sig_synth;
    
    carrier_synth(:,i) = sig_synth;
end
%clear sinu res F M P mind F_ sinf M_ sinm sin_synth coeff ar_model res_synth sig_synth
toc

demod_struct.carriers = carriers;
demod_struct.carrier_synth = carrier_synth;

end