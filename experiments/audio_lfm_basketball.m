addpath(genpath('../lfm_augmented'))
addpath(genpath('../lfm_toolbox'))
addpath(genpath('../sinemodel'))
addpath(genpath('../gppad'))
addpath(genpath('../filterbank'))


seed = 12345;
rng(seed,'twister') % for repeatability

fileName = 'stim18_basketball_dribbling'; % name of wav file
suffix = ''; % suffix for .mat file to be saved (in case of multiple versions)
saveResults = 1;    % save the results to a mat file?
matNameLFM = strcat('../audio_lfm_data/lfm_',fileName,suffix,'.mat'); % name of the .mat file saved if saveResults = 1
matNameDemod = strcat('../audio_lfm_data/demod_',fileName,suffix,'.mat'); % name of demodulation .mat file
numFilt = 16; % number of subbands
len = [100 100]; % range of demodulation lengthscales (different frequencies tend to need different lengthscales)
loadPrevDemod = 1; % load previously saved demodulation results
loadPrevLFM = 1; % load previously saved latent force model results
run_lfm = 0; % run the latent force model
run_append = 0; % append extra envelopes to the initial LFM
run_final = 0; % run the final filtering stage with all envelopes
run_post = 0; % run post-processing of latent force
hil4lo = 0; % use the Hilbert envelope for this many subbands (starting from lowest)

if loadPrevDemod == 0
    demod_struct = demod_filterbank(fileName,suffix,numFilt,loadPrevDemod,len,hil4lo);
    save(matNameDemod,'-struct','demod_struct');
else
    demod_struct = load(matNameDemod);
end

envs = demod_struct.A'; % envelopes / modulators

analysisRange = 1:12400; % analysis range (chosen manually - smallest representative chunk of the signal)
DS = 50; % downsample rate
N = 6; % # envelopes to include in first (main) LFM (all other envelopes are appended to the model optimised on the N envelopes with the greatest energy)

lfmOpt = struct; % store the LFM settings in a structure
lfmOpt.analysisRange = analysisRange;
lfmOpt.DS = DS;
lfmOpt.nlf = 1; % # latent forces
lfmOpt.sfactor = 1; % # lfm downsample rate (not implemented)
lfmOpt.priorityInference = 0; % give priority to high amp. outputs (not implemented)
lfmOpt.dt = 0.1; % time step size
lfmOpt.gamma = 0.99; % linearity measure
lfmOpt.grng = [0.8 1]; % allowable range for gamma
lfmOpt.Dinit = 1.25; % initial guess for damping
lfmOpt.Drng = [0.75 4.5]; % allowable range of damping values
lfmOpt.linit = 1; % initial guess for lengthscale
lfmOpt.lrng = [0.12 20]; % allowable range of lengthscales
lfmOpt.Srng = [0 0.4]; % allowable range for sensitivity
lfmOpt.fbinit = 0.02; % initial guess for feedback terms
lfmOpt.fwdinit = 0.02; % initial guess for forward terms
lfmOpt.sig_noise = 0.00175; % assumed measurement noise variance (automatically calculated below)
lfmOpt.silenceThreshold = 0.05; % if maxmimum amplitude drops below this value times the max of the signal, don't make predictions
lfmOpt.silThrCovReset = 0; % if maxmimum amplitude drops below this value times the max of the signal, rest the covariance to the prior
lfmOpt.resetOpt = 0; % perform covariance resetting in silent periods during optimisation?
lfmOpt.resetFinal = 0; % perform covariance resetting in silent periods during final filtering and smoothing?
lfmOpt.preset = 0; % periodic reset time
lfmOpt.preset_it = 0; % number of iterations to perform periodic reset for
%lfmOpt.restrictFunEvals = 1; % 1 if we restrict number of function evaluations
    lfmOpt.MaxFunEvals = 6000;
%lfmOpt.restrictIter = 1; % 1 if we restrict number of iterations
    lfmOpt.MaxIter = 100;
    lfmOpt.TolFun = 0.0015;
    
% FEEDBACK
lfmOpt.fb_terms = 18; % number of additional feedback terms
lfmOpt.sparsity_fb = [0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 1]; % index which feedback terms to include (omit to optimise all)
lfmOpt.lambda_damping = 0; % regulariser for damping magnitudes (L2 norm) - helps ensure stability
lfmOpt.lambda_sensitivity = 0; % regulariser for sensitivity (L2 norm) - helps constrain / localise optimisation
lfmOpt.lambda_feedback = 0; % regulariser for feedback terms (L1 norm) - required if optimisation is not sparse
%lfmOpt.threshold_up_damping = 5;%3.5; % punish damping coefficients over this (soft) threshold
%lfmOpt.threshold_lo_damping = 0;%0.75; % punish damping coefficients under this (soft) threshold
%lfmOpt.threshold_lengthscale = 0;%0.12; % punish lengthscales under this (soft) threshold

% FORWARD
lfmOpt.fwd_terms = 24; % number of additional forward terms
lfmOpt.sparsity_fwd = [0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 1]; % index which feedback terms to include (omit to optimise all)
lfmOpt.lambda_lag = 0; % regulariser for distributed lag (L1 norm) - required if optimisation is not sparse

lfmData = struct;
lfmData.fileName = fileName;
lfmData.suffix = suffix;
lfmData.demod_struct = demod_struct;
lfmData.lfmOpt = lfmOpt;
lfmData.N = N;
lfmData.modulators = envs;


% Sort subbands by energy contribution
maxEn = max(envs,[],2) / max(max(envs,[],2));
totEn = sum(envs,2) / max(sum(envs,2));
En = maxEn + totEn;
[~,I] = sort(En,'descend');
lfmData.Energy = En;
lfmData.I = I;
envs_sorted = envs(I,:);
% adjust so max value is 0.25:
env_scale = 0.25 / max(max(envs));
envelopes = envs_sorted * env_scale; % the envelopes to be analysed by the LFM

lfmData.envs_sorted = envs_sorted;
lfmData.env_scale = env_scale;
lfmData.envelopes = envelopes;

% plot the modulators
figure(2); clf
subplot(411)
plot(envs')
title('The Amplitude Envelopes')
subplot(412)
plot(analysisRange,envs(:,analysisRange)')
title('A Small Range of the Envelopes')
subplot(413)
plot(envs(:,analysisRange(1):DS:analysisRange(end))')
title(sprintf('Downsampled by factor of %d',DS))
subplot(414)
plot(envelopes(1:N,analysisRange(1):DS:analysisRange(end))')
title('Envelopes with Greatest Energy (Scaled)')
drawnow();

sprev = rng(); % for repeatability

%% Space for manual editing of optimisation parameters:


%% load previously run version
theta_guess = [];
if loadPrevLFM == 1
    load(matNameLFM);
    theta_guess = lfmData.lfmResults.theta_opt;
end
%% Run the latent force model
if run_lfm == 1
    rng(sprev)
    fprintf('Running %s: %d outputs with %d latent force(s)\n',strcat(fileName,suffix),N,lfmData.lfmOpt.nlf)
    % train on the selected range (slow)
    results = runLFM(lfmData,0,theta_guess,1,[],1:N);
    lfmData.lfmResults = results;
    if saveResults == 1
        save(matNameLFM,'lfmData');
    end
    % then filter and smooth entire signal with learnt parameters (faster)
    filterSig = runLFM(lfmData,1,lfmData.lfmResults.theta_opt,0,[],1:N);
    lfmData.filterSig = filterSig;
    % save
    if saveResults == 1
        save(matNameLFM,'lfmData');
    end
    drawnow();
end

%% Now fix the learnt parameters and append the smaller envelopes
if run_append == 1
    theta_guess = lfmData.lfmResults.theta_opt;
    %fixDim = 3:4 % TODO: ADD FUNCTIONALITY TO FIX ARBITRARY ENVELOPES (NOT JUST THE FIRST N)
    appendN = N-2;
    lfmData.appendN = appendN;
    fixDim = 1:appendN;
    numNewEnvs = 2; % # envelopes to append at a time
    optDim = N;
    append_results = cell(ceil((size(envs,1)-N)/numNewEnvs),1);
    lfmData.lfmOpt.MaxIter = 15;
    lfmData.lfmOpt.TolFun = 0.002;
    for i=1:ceil((size(envs,1)-N)/numNewEnvs)
        optDim = optDim(end)+1:min(optDim(end)+numNewEnvs,size(envs,1));
        fprintf('Running %s: Appending %d output(s): %s\n',strcat(fileName,suffix),length(optDim),num2str(optDim))
        append_results{i} = runLFM(lfmData,0,theta_guess,1,fixDim,optDim);
        append_results{i}.fixDim = fixDim;
        append_results{i}.optDim = optDim;
        lfmData.append_results = append_results;
        if saveResults == 1
            save(matNameLFM,'lfmData');
        end
    end
end

%% Combine results into one parameter vector and run the final filtering stage with all the envelopes
if run_final == 1
    [ls, gamma, D, ~, ~, S_sens, S_fb_, S_fwd_] = unpackTheta(lfmData.filterSig);
    S_fb = S_fb_(:,logical(lfmData.lfmOpt.sparsity_fb));
    sparsity_fwd_ind = [];
    for i=1:lfmData.lfmOpt.nlf
        sparsity_fwd_ind = [sparsity_fwd_ind lfmData.lfmOpt.sparsity_fwd]; %#ok<AGROW>
    end
    S_fwd = S_fwd_(:,logical(sparsity_fwd_ind));

    for i=1:size(lfmData.append_results,1)
        [~, gamma_append, D_append, ~, ~, S_sens_append, S_fb_append, S_fwd_append] = unpackTheta(lfmData.append_results{i});
        S_fb_append = S_fb_append(:,logical(lfmData.lfmOpt.sparsity_fb));
        S_fwd_append = S_fwd_append(:,logical(sparsity_fwd_ind));
        gamma = [gamma; gamma_append(lfmData.appendN+1:end)]; %#ok<AGROW>
        D = [D; D_append(lfmData.appendN+1:end)]; %#ok<AGROW>
        S_sens = [S_sens; S_sens_append(lfmData.appendN+1:end,:)]; %#ok<AGROW>
        S_fb = [S_fb; S_fb_append(lfmData.appendN+1:end,:)]; %#ok<AGROW>
        S_fwd = [S_fwd; S_fwd_append(lfmData.appendN+1:end,:)]; %#ok<AGROW>
    end
    theta_final = [ls; gamma; D; S_sens(:); S_fb(:); S_fwd(:)];
    final_subset = 1:6;%numFilt;
    theta_final_subset = [ls; gamma(final_subset); D(final_subset); reshape(S_sens(final_subset,:),size(S_sens,2)*length(final_subset),1); ...
                                reshape(S_fb(final_subset,:),size(S_fb,2)*length(final_subset),1); ...
                                    reshape(S_fwd(final_subset,:),size(S_fwd,2)*length(final_subset),1)];
    lfmData.theta_final = theta_final;
    lfmData.theta_final_subset = theta_final_subset;
    % then filter and smooth entire signal with all the envelopes
    fprintf('Running %s: final filtering stage with all %d outputs and %d latent force(s)\n',strcat(fileName,suffix),size(lfmData.envelopes,1),lfmData.lfmOpt.nlf)
    filterFinal = runLFM(lfmData,1,theta_final_subset,0,[],final_subset);
    lfmData.filterFinal = filterFinal;
    
    lfmData.filterFinal.theta_opt = theta_final;
    lfmData.filterFinal.inSig = lfmData.envelopes(:,1:lfmData.lfmOpt.DS:end);
    lfmData.filterFinal.At = lfmData.filterFinal.inSig(:,1);
    lfmData.filterFinal.Bt = zeros(size(lfmData.filterFinal.At));
    
    % save
    if saveResults == 1
        save(matNameLFM,'lfmData');
    end
end

%% Post-processing of the force to deal with potential noise floor and improve Euclidean distance
if run_post == 1
    fprintf('Post processing of the latent force\n'); tic
    lfmData = forceProcessing(lfmData,5); toc
    if saveResults == 1
        save(matNameLFM,'lfmData');
    end
end

%% Resynthesis

RemoveNoiseFloor = 0; % set to >0 to remove noise floor
playSounds = 0;

[envSynth,latentForces,sigSynth] = reconstruct_sig(lfmData,[],RemoveNoiseFloor);

lfmData.envSynth = envSynth;
lfmData.latentForces = latentForces;
lfmData.sigSynth = sigSynth;

actEnv = lfmData.demod_struct.A';
[w,h] = nnmf(actEnv,lfmData.lfmOpt.nlf); % NMF
nmfSynth = w*h;

[tnmf_w,tnmf_h] = run_tNMF(actEnv(:,1:lfmData.lfmOpt.DS:end),lfmData.lfmOpt.nlf); % temporal NMF
tnmf_h = makeThisThatLong(tnmf_h,size(actEnv,2)); tnmf_h(find(tnmf_h<0))=0; %#ok<FNDSB>

tnmfSynth = tnmf_w*tnmf_h;
tnmfSynth_a = bsxfun(@times,tnmfSynth,max(actEnv,[],2)./max(tnmfSynth,[],2)); % magnitude adjustment
if sum(sum((tnmfSynth_a-actEnv).^2)) < sum(sum((tnmfSynth-actEnv).^2));
    tnmfSynth = tnmfSynth_a;
end

errLFM = sum(sum((envSynth'-actEnv).^2));
errNMF = sum(sum((nmfSynth-actEnv).^2));
errTNMF = sum(sum((tnmfSynth-actEnv).^2));
fprintf('LFM error: %s\n',errLFM)
fprintf('NMF error: %s\n',errNMF)
fprintf('tNMF error: %s\n',errTNMF)

lfmdist = 0;
nmfdist = 0;
tnmfdist = 0;
for i=1:size(lfmData.envelopes,1)
    lfmdist = lfmdist + pdist([envSynth(:,i)';actEnv(i,:)],'cosine');
    nmfdist = nmfdist + pdist([nmfSynth(i,:);actEnv(i,:)],'cosine');
    tnmfdist = tnmfdist + pdist([tnmfSynth(i,:);actEnv(i,:)],'cosine');
end
fprintf('LFM cos dist.: %s\n',lfmdist)
fprintf('NMF cos dist.: %s\n',nmfdist)
fprintf('tNMF cos dist.: %s\n',tnmfdist)

lfmData.objective_results.errLFM = errLFM;
lfmData.objective_results.errNMF = errNMF;
lfmData.objective_results.errTNMF = errTNMF;
lfmData.objective_results.lfmdist = lfmdist;
lfmData.objective_results.nmfdist = nmfdist;
lfmData.objective_results.tnmfdist = tnmfdist;
if saveResults == 1
    save(matNameLFM,'lfmData');
end

nmfSig = sum(nmfSynth' .* demod_struct.carrier_synth,2);

if playSounds == 1
    sound(lfmData.demod_struct.y,lfmData.demod_struct.fs)
    pause(length(lfmData.demod_struct.y)/lfmData.demod_struct.fs + 0.5)
    sound(sigSynth,lfmData.demod_struct.fs)
end

figure(6);clf
subplot(331)
plot(actEnv')
title(strcat('Actual',{' '},fileName,{' '},'Envelopes'),'interpreter','none')
yLimits = ylim;
subplot(334)
plot(latentForces','LineWidth',2)
title('LFM Latent Function(s)')
subplot(337)
plot(envSynth)
ylim(yLimits)
title('LFM Reproduction')
annotation('textbox',[0.235 0.21 0.1 0.1],'FontSize',12,'FontWeight','Bold','String',strcat('Euclidean Dist.:',{' '},num2str(errLFM)),'FitBoxToText','on', 'EdgeColor', 'none','Color','k')
annotation('textbox',[0.235 0.185 0.1 0.1],'FontSize',12,'FontWeight','Bold','String',strcat('Cosine Dist.:',{' '},num2str(lfmdist)),'FitBoxToText','on', 'EdgeColor', 'none','Color','k')
subplot(332)
plot(actEnv')
title(strcat('Actual',{' '},fileName,{' '},'Envelopes'),'interpreter','none')
subplot(335)
plot(h','LineWidth',2)
title('NMF Latent Function(s)')
subplot(338)
plot(nmfSynth')
ylim(yLimits)
title('NMF Reproduction')
annotation('textbox',[0.515 0.21 0.1 0.1],'FontSize',12,'FontWeight','Bold','String',strcat('Euclidean Dist.:',{' '},num2str(errNMF)),'FitBoxToText','on', 'EdgeColor', 'none','Color','k')
annotation('textbox',[0.515 0.185 0.1 0.1],'FontSize',12,'FontWeight','Bold','String',strcat('Cosine Dist.:',{' '},num2str(nmfdist)),'FitBoxToText','on', 'EdgeColor', 'none','Color','k')
subplot(333)
plot(actEnv')
title(strcat('Actual',{' '},fileName,{' '},'Envelopes'),'interpreter','none')
subplot(336)
plot(tnmf_h','LineWidth',2)
title('tNMF Latent Function(s)')
subplot(339)
plot(tnmfSynth')
ylim(yLimits)
title('tNMF Reproduction')
annotation('textbox',[0.795 0.21 0.1 0.1],'FontSize',12,'FontWeight','Bold','String',strcat('Euclidean Dist.:',{' '},num2str(errTNMF)),'FitBoxToText','on', 'EdgeColor', 'none','Color','k')
annotation('textbox',[0.795 0.185 0.1 0.1],'FontSize',12,'FontWeight','Bold','String',strcat('Cosine Dist.:',{' '},num2str(tnmfdist)),'FitBoxToText','on', 'EdgeColor', 'none','Color','k')