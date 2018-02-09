function genAnchor(f_ind, sig_dur, suffix)

% function to generate the anchor (i.e. poor quality reference sound) for
% the MUSHRA-style listening test

% Will Wilkinson Jan-2018

if nargin < 3
    suffix = '';
end

addpath(genpath('GpMat'))
addpath(genpath('../gppad'))

% List of all examples
fileNames = {'metal','bike_horn2','stim360_cymbal_crash','stim186_hammering','gunshot2','squeaky_stairs2','squeaky_stairs3',... % 1-7
             'stim7_applause','stim41_breathing','stim162_frog_croaking','glass','dog','stim117_dog_drinking','stim114_dishes_clanking',... % 8-14
             'stim18_basketball_dribbling','stim373_train_warning_bell','stim132_water_dripping','stim373_train_warning_bell2','rain','stim399_walking_on_hard_surface','stim78_chimes_in_the_wind',... % 15-21
             'stim211_keys_jingling','stim312_wind','stim50_camera_snapping_photos','cello','stim398_walking_on_gravel'}; % 22 - end

%f_ind = 1; % file name index, e.g. 8 = applause
fileName = fileNames{f_ind};

fprintf('generating anchor for %s\n',fileName);

% Load LFM results into a struct called 'lfmData'
load(strcat('../audio_lfm_data/lfm_',fileName,'.mat'));

nlf = 1;%lfmData.lfmOpt.nlf;
dt = lfmData.lfmOpt.dt;

%sig_dur = 4; % duration of synthetic signal in seconds
dur_in_samp = floor(sig_dur * lfmData.demod_struct.fs / lfmData.lfmOpt.DS);
playSound = 1; % play the synthetic sound?
saveSound = 1; % save the synthetic sound?
useLFM = 0; % 0 = NMF, 1 = LFM

%%
%seed = 0987;
%rng(seed,'twister') % for repeatability

ds = lfmData.lfmOpt.DS;
[nmf_w,nmf_h] = nnmf(lfmData.demod_struct.A(1:ds:end,:)',nlf); % NMF

nmf_h(find(nmf_h<=0.0001)) = 0.0001; %#ok<FNDSB>

%%

h_f = @(f) log(exp(f)-1);
g_f = @(f) log(1+exp(f));

if useLFM
    forceC = log(1+exp(lfmData.filterFinal.uAdjusted));
else
    forceC = nmf_h;
end

forceC_neg = h_f(forceC);
fCmu = mean(forceC_neg,2);
forceC_adjust = bsxfun(@minus,forceC_neg,fCmu);

for n=1:nlf
    forceC_adjust(n,find(forceC_adjust(n,:)<-max(forceC_adjust(n,:)))) = -max(forceC_adjust(n,:)); %#ok<FNDSB>
end
forceC_adjust = bsxfun(@minus,forceC_adjust,mean(forceC_adjust,2));

figure(1);clf
subplot(211)
plot(nmf_h')
subplot(212)
plot(forceC_adjust')
drawnow();

%%
DS = 2;
% Discretize x (measurement points)
%x = linspace(1,size(fA_tr,1),ceil(size(fA_tr,1)/DS))';
% Test points (evaluation points)
xt = linspace(1,size(nmf_h,2),size(nmf_h,2))';
% Measurement noise variance
sigma2 = .001;

%% Train GP on latent forces to get better parameters for the generative model
fprintf('optimising latent function parameters\n');
xfC_ind = cell(nlf,1);
xfC_logicals = bsxfun(@gt,forceC_adjust',-prctile(forceC_adjust',90));
newLF = cell(nlf,1);
forceC_neg_train = forceC_neg;
for fnum = 1:nlf
    
    xfC = find(xfC_logicals(:,fnum));
    xfC_ind{fnum} = xfC;
    x = xfC_ind{fnum}(1:DS:end);
    
    forceC_neg_train(fnum,xfC_ind{fnum}) = NaN;
    
    obsv = forceC_adjust(fnum,:)'; % observations

      % Simulate data
      y = obsv(x) + sqrt(sigma2)*randn(size(x));

    % Perform inference
      % Noise variance prior
      ps2 = prior_logunif(); 

      % The likelihood model
      lik = lik_gaussian('sigma2', 1, 'sigma2_prior', ps2);

      % Covariance function hyperparameter priors
      pl = prior_logunif(); 
      pm = prior_logunif();

      % The GP covariance function
      gpcf = gpcf_sexp('lengthScale', 1, ...
                           'magnSigma2', 1, ...
                           'lengthScale_prior', pl, ...
                           'magnSigma2_prior', pm);

      % Define Gaussian process model using type 'KALMAN'
      gp = gp_set('lik', lik, 'cf', gpcf, 'type', 'KALMAN');

      % Hyperparameter optimization (state space model)
      warning('off','all');
      gp = gp_optim(gp, x, y);
      warning('on','all');
      
      newLF{fnum}.lengthScale = gp.cf{1}.lengthScale;
      newLF{fnum}.magnSigma2 = gp.cf{1}.magnSigma2;
      newLF{fnum}.u = obsv;
      newLF{fnum}.mu = fCmu(fnum);
end
    
%% Sample from the prior and plot

figure(3); clf
subplot(221)
plot(forceC_neg')
hold on
plot(forceC_neg_train','k')
xlim1= xlim; ylim1 = ylim;
subplot(223)
plot(forceC')
xlim2 = xlim; ylim2 = ylim;
lowLevelGen = zeros(dur_in_samp,nlf);
scaleSigma = ones(1,nlf);%[1.75 4.25 3.5];%ones(1,nlf);
sparseParam = -0*ones(1,nlf);%[-1 -1 -2];%-0*ones(1,nlf);
for fnum=1:nlf
    fprintf('sampling latent function %d\n',fnum);
    warning('off','all');
    % multiply lengthscale by 2 to make the anchor intentionally too
    % smooth
    fsamp = gpSample('rbf', 1, [2*1/newLF{fnum}.lengthScale^2 scaleSigma(fnum)*newLF{fnum}.magnSigma2], [1 dur_in_samp]);%size(newLF{fnum}.u,1)]);
    warning('on','all');
    fsamp_tr = g_f(fsamp + newLF{fnum}.mu + sparseParam(fnum));
    fsamp_tr = (max(nmf_h')/max(fsamp_tr))*fsamp_tr;
    % shift down to make the anchor intentionally too sparse
    fsamp_tr = fsamp_tr - 0.05*max(fsamp_tr);
    fsamp_tr(find(fsamp_tr<=0)) = 0.00001; %#ok<FNDSB>
    newLF{fnum}.fsamp = fsamp;
    newLF{fnum}.fsamp_tr = fsamp_tr;
    lowLevelGen(:,fnum) = fsamp_tr;
    subplot(222)
    plot(fsamp + newLF{fnum}.mu)
    hold on
    axis([1 Inf ylim1])
    subplot(224)
    plot(fsamp_tr)
    hold on
    axis([1 Inf ylim2])
end
drawnow();

% Multiply cascade levels to produce envelopes
latentGen = lowLevelGen;
figure(4);clf
subplot(221)
plot(nmf_h')
xlim_in = xlim; ylim_in = ylim;
subplot(223)
plot(latentGen)
axis([1 Inf ylim_in])
drawnow();

if useLFM
    [envSynth,~,sigSynth] = reconstruct_sig(lfmData,latentGen',[],[],sig_dur*lfmData.demod_struct.fs,[],[],zeros(lfmData.demod_struct.D,1));
else
    envSynth = nmf_w * latentGen';
    envSynth = makeThisThatLong(envSynth,sig_dur*lfmData.demod_struct.fs)';

    carrier_synth = lfmData.demod_struct.carrier_synth;
    carriers = zeros(size(envSynth));
    for i=1:size(envSynth,1)
        carriers(i,:) = carrier_synth(mod(i-1,size(carrier_synth,1))+1,:);
    end
    
    sigSynth = sum(carriers'.*envSynth');
end

figure(4)
subplot(222)
plot(lfmData.demod_struct.A)
title('Actual Envelopes')
xlim_out = xlim; ylim_out = ylim;
subplot(224)
plot(envSynth)
title('Reproduced Envelopes')
axis([1 Inf ylim_out])
drawnow();

if playSound == 1
    sound(sigSynth,lfmData.demod_struct.fs)
end
if saveSound == 1
    audiowrite(strcat('../audio/synth/',fileName,'_anchor',suffix,'.wav'),sigSynth,lfmData.demod_struct.fs);
end

end