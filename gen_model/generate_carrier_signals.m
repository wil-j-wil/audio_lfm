addpath(genpath('../sinemodel'))
addpath(genpath('../filterbank'))

% List of all examples
fileNames = {'metal','bike_horn2','stim360_cymbal_crash','stim186_hammering','gunshot2','squeaky_stairs2','squeaky_stairs3',... % 1-7
             'stim7_applause','stim41_breathing','stim162_frog_croaking','glass','dog','stim117_dog_drinking','stim114_dishes_clanking',... % 8-14
             'stim18_basketball_dribbling','stim373_train_warning_bell','stim132_water_dripping','stim373_train_warning_bell2','rain','stim399_walking_on_hard_surface','stim78_chimes_in_the_wind',... % 15-21
             'stim211_keys_jingling','stim312_wind','stim50_camera_snapping_photos','cello','stim398_walking_on_gravel'}; % 22 - end

%f_ind = 1;
for f_ind = 1:26
    
    fileName = fileNames{f_ind};
    % Load LFM results into a struct called 'lfmData'
    load(strcat('../audio_lfm_data/lfm_',fileName,'.mat'));

    demod_struct = lfmData.demod_struct;
    C = demod_struct.C;
    D = demod_struct.D;
    Fc = demod_struct.Fc;
    y = demod_struct.y;
    fs = demod_struct.fs;

    sig_dur = 4; % in seconds
    dur_in_samp = ceil(sig_dur * fs);

    %% Analyse Subband Carriers
    carriers = cell(D,1);
    AR_order = 120;
    envLen = max(350-floor(Fc),50); % length of smoothing/normalisation envelope
    carrier_synth = zeros(dur_in_samp,size(C,2));
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

        sinf = zeros(size(M,1),1);
        sinm = zeros(size(M,1),1);
        sin_synth = zeros(size(M,1),dur_in_samp);
        for j=1:size(M,1)
            F_ = F(j,:);
            M_ = M(j,:);
            sinf(j) = median(F_(find(isnan(F_)==0))); %#ok<FNDSB>
            sinm(j) = median(M_(find(isnan(M_)==0))); %#ok<FNDSB>
            sin_synth(j,:) = sinm(j)*sin(2*pi*sinf(j)*((1:dur_in_samp)+rand(1))/fs); % random initial phase
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

        coeff = lpc(res,AR_order);
        ar_model = arima('Constant',0,'AR',-coeff(2:end),'Variance',.00001);
        res_synth = simulate(ar_model,dur_in_samp);

        res_synth = median(envelope(res,envLen(i),'rms')) * res_synth ./ median(envelope(res_synth,envLen(i),'rms'));

        if isempty(sin_synth)
            sig_synth = res_synth;
        else
            sig_synth = sin_synth + res_synth;
        end

        carrier_synth(:,i) = sig_synth;
    end
    toc

    %% Save
    if ~isfield(lfmData.demod_struct,'carrier_synth_old')
        lfmData.demod_struct.carrier_synth_old = lfmData.demod_struct.carrier_synth;
    end
    lfmData.demod_struct.carrier_synth = carrier_synth;
    save(strcat('../audio_lfm_data/lfm_',fileName,'.mat'),'lfmData');
    
end