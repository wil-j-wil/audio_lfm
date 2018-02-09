function [cascadeMod, newLF] = genModelTNMF(f_ind, sig_dur, suffix)

% function to generate novel sounds using the tNMF-based model

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
         
lengthscales = {   5,      5,   10,   20,   10,   20,   13,...
                16.6,     20,   50,    5,   10,   30,   11,...
                  20,      4,   25,    5,   75,  100, 11.6,...
                  30,   31.2,   50, 16.6,   50};
            
fA_prctile   = {60, 20, 60, 30, 60, 20, 20,...
                20, 20, 20, 20, 20, 20, 20,...
                20, 20, 00, 20, 20, 40, 20,...
                20, 40, 00, 60, 20};
            
scale_sparse = { [1.2 3], [0.5 10],   [10 3],     [1 7],   [1.2 4],      [2 5],  [1.4 6],...
                [0.8 50], [0.9 10],    [1 5],  [1.5 10],     [4 5],  [0.95 10], [1.5 10],...
                   [3 2],  [1 100], [4.75 3],   [2 2.5],   [1 100], [0.95 100],    [5 3],...
                 [1.2 5],  [1 100],  [1 100],  [0.5 10], [1.25 10]};

%f_ind = 26; % file name index, e.g. 8 = applause
fileName = fileNames{f_ind};

fprintf('generating tNMF sample for %s\n',fileName);

% Load LFM results into a struct called 'lfmData'
load(strcat('../audio_lfm_data/lfm_',fileName,'.mat'));

nlf = lfmData.lfmOpt.nlf;
dt = lfmData.lfmOpt.dt;
fs = lfmData.demod_struct.fs;

%sig_dur = 4; % duration of synthetic signal in seconds
dur_in_samp = floor(sig_dur * fs / lfmData.lfmOpt.DS);
playSound = 0; % play the synthetic sound?
saveSound = 1; % save the synthetic sound?

[T,D] = size(lfmData.demod_struct.A);
%T=20000;
ATrain = lfmData.demod_struct.A(1:lfmData.lfmOpt.DS:T,:);
T = size(ATrain,1);

%% Train the model
K = nlf;
HInit = exp(randn(T,K));
ks = ceil(T*rand(K,1));
WInit = ATrain(ks,:);
vary = zeros(T,D);

% initialise with regular NMF
Opts.restarts = 50;
Opts.numIts = 5000;
%[WEst1,HEst1,info1] = nmf(ATrain,WInit,HInit,[],[],[],vary,Opts);
[WEst1,HEst1,info1] = nmf_fp(ATrain,WInit,HInit,vary,Opts);

% order according to slowness (mean square derivative)
fastness = mean(diff(HEst1).^2)./var(HEst1);
[val,ind] = sort(fastness,'descend');
HEst1 = HEst1(:,ind); 
WEst1 = WEst1(ind,:);

%% LOG-NORMAL SE GAUSSIAN PROCESS TEMPORAL NMF
% initialise

lenx = zeros(K,1);
mux = zeros(K,1);
varx = zeros(K,1);

for k=1:K
  % threshold
  logHthresh = log(HEst1(:,k)+1e-8);
  
  % filter
  filt = exp(-1/2*([-100:100].^2)/(1^2));
  filt = filt/sum(filt);
  logHsm = conv(logHthresh,filt,'same');
  
  % fit GP
  mux(k) = mean(logHsm);
  [lenx(k),varx(k),info] = trainSEGP_RS(logHsm-mux(k));
  
end

disp(['length scale parameters'])
lenx
Opts.numIts = 200;
[HEst2,info2] = tnmf_inf(ATrain,WEst1,HEst1+1e-8,lenx,mux,varx,vary,Opts);
Opts.numIts = 2000;
[WEst2,HEst2,info2] = tnmf(ATrain,WEst1,HEst2,lenx,mux,varx,vary,Opts);

%%
force_pos = HEst2;

if nlf>1
    force_max = max(force_pos');
else
    force_max = force_pos';
end
%len = lfmData.demod_struct.len(end)/50;
%len = max(2.5,2.5*max(lfmData.lfmResults.theta_opt(1:nlf))/lfmData.lfmOpt.dt);
len = lengthscales{f_ind};

h_f = @(f) log(exp(f)-1);
g_f = @(f) log(1+exp(f));

[fA,fC] = GPPAD(force_max',len,0);
forceC = bsxfun(@rdivide,HEst2,fA)';

fA_neg = h_f(fA);
fAmu = mean(fA_neg);
fA_tr = fA_neg - fAmu;

forceC_neg = h_f(forceC);
fCmu = mean(forceC_neg,2);
forceC_adjust = bsxfun(@minus,forceC_neg,fCmu);

for n=1:nlf
    forceC_adjust(n,find(forceC_adjust(n,:)<-max(forceC_adjust(n,:)))) = -max(forceC_adjust(n,:)); %#ok<FNDSB>
end
forceC_adjust = bsxfun(@minus,forceC_adjust,mean(forceC_adjust,2));

figure(1);clf
subplot(511)
plot(HEst2)
subplot(512)
plot(force_max','b')
hold on
plot(fA,'r')
subplot(513)
plot(fA_tr)
subplot(514)
plot(forceC_adjust')
subplot(515)
plot(forceC')
drawnow();

%%
DS = 2;
% Discretize x (measurement points)
%x = linspace(1,size(fA_tr,1),ceil(size(fA_tr,1)/DS))';
% Test points (evaluation points)
xt = linspace(1,size(fA_tr,1),size(fA_tr,1))';
% Measurement noise variance
sigma2 = .001;


%% Train GP on high-level cascade modulator
fprintf('optimising high-level modulator parameters\n');
% only train on regions of high energy (upper 2 quartiles)
%x = find(fA>prctile(fA,50));
xfA_ind = find(fA_tr>max(-prctile(fA_tr,99),prctile(fA_tr,fA_prctile{f_ind})));
x = xfA_ind(1:DS:end);
%x = x(1:DS:end);
%x = linspace(1,size(fA_tr,1),ceil(size(fA_tr,1)/DS))';

obsv = fA_tr; % observations

% Simulate data
%y = obsv(1:DS:end) + sqrt(sigma2)*randn(size(x));
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

cascadeMod.lengthScale = gp.cf{1}.lengthScale;
cascadeMod.magnSigma2 = gp.cf{1}.magnSigma2;
cascadeMod.u = obsv;
cascadeMod.mu = fAmu;

%% Train GP on latent forces to get better parameters for the generative model
fprintf('optimising latent function parameters\n');
xfC_ind = cell(nlf,1);
xfC_logicals = bsxfun(@gt,forceC_adjust',-prctile(forceC_adjust',90));
newLF = cell(nlf,1);
forceC_neg_train = forceC_neg;
for fnum = 1:nlf
    
    xfC = find(xfC_logicals(:,fnum));
    xfC_ind{fnum} = intersect(xfC,xfA_ind);
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
for suf = 1:length(suffix)
    
    fA_tr_train = fA_tr;
    fA_tr_train(xfA_ind) = NaN;

    figure(2); clf
    subplot(411)
    plot(fA)
    xlim0= xlim; ylim0 = ylim;
    subplot(412)
    plot(fA_tr)
    hold on
    plot(fA_tr_train,'k')
    
    fprintf('calculating sparsity parameters\n');
    warning('off','all');
    fsamp_temp = gpSample('rbf', 5, [1/cascadeMod.lengthScale^2 cascadeMod.magnSigma2], [1 2500]);%size(cascadeMod.u,1)]);
    warning('on','all');

    u_pos = g_f(cascadeMod.u);
    u_normalized = u_pos/max(abs(u_pos));
    eU = entropy(u_normalized);
    gp_pos = g_f(fsamp_temp);
    gp_normalized = bsxfun(@rdivide,gp_pos,max(abs(gp_pos),[],2));
    eGP = entropy(gp_normalized);

    sparseParam = min(max(eGP - eU,-0.2),2)
    scaleSigma = 1 + sparseParam / cascadeMod.magnSigma2

    %scaleSigma = 2.6;%scale_sparse{f_ind}(1);%1-cascadeMod.mu/20;%1
    %sparseParam = -2;%-abs(cascadeMod.mu)/scale_sparse{f_ind}(2);%-0
    fprintf('sampling high-level modulator\n');
    warning('off','all');
    fsamp = gpSample('rbf', 1, [1/cascadeMod.lengthScale^2 scaleSigma*cascadeMod.magnSigma2], [1 dur_in_samp]);%size(cascadeMod.u,1)]);
    fsamp_tr = g_f(fsamp + cascadeMod.mu - sparseParam);
    if max(fsamp_tr)>1.6*max(fA) || max(fsamp_tr)<0.75*max(fA)
        fprintf('resampling\n')
        fsamp = gpSample('rbf', 1, [1/cascadeMod.lengthScale^2 scaleSigma*cascadeMod.magnSigma2], [1 dur_in_samp]);%size(cascadeMod.u,1)]);
        fsamp_tr = g_f(fsamp + cascadeMod.mu - sparseParam);
    end
    warning('on','all');
    %fsamp_pre = fsamp;
    %if max(fsamp)>1.5*max(fA_tr)
    %    fsamp = (1.5*max(fA_tr)/max(fsamp))*fsamp;
    %elseif max(fsamp)<0.9*max(fA_tr)
    %    fsamp = (0.9*max(fA_tr)/max(fsamp))*fsamp;
    %end
    %fsamp_tr = g_f(fsamp + cascadeMod.mu - sparseParam);
    if max(fsamp_tr)>1*max(fA)
        fsamp_tr = (1*max(fA)/max(fsamp_tr))*fsamp_tr;
    elseif max(fsamp_tr)<0.5*max(fA)
        fsamp_tr = (0.5*max(fA)/max(fsamp_tr))*fsamp_tr;
    end
    cascadeMod.fsamp = fsamp;
    cascadeMod.fsamp_tr = fsamp_tr;
    subplot(413)
    plot(fsamp)
    subplot(414)
    hold on
    plot(g_f(fsamp + cascadeMod.mu - sparseParam),'r--')
    plot(fsamp_tr)
    axis([1 Inf ylim0])
    drawnow();

    maxForce = prctile(forceC',99)';

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
    %scaleSigma = ones(1,nlf);%[1.75 4.25 3.5];%ones(1,nlf);
    %sparseParam = -0*ones(1,nlf);%[-1 -1 -2];%-0*ones(1,nlf);
    for fnum=1:nlf
        fprintf('sampling latent function %d\n',fnum);
        scaleSigma = 0.92;%1-newLF{fnum}.mu/20;
        sparseParam = newLF{fnum}.mu/10;
        warning('off','all');
        fsamp = gpSample('rbf', 1, [1/newLF{fnum}.lengthScale^2 scaleSigma*newLF{fnum}.magnSigma2], [1 dur_in_samp]);%size(newLF{fnum}.u,1)]);
        warning('on','all');
        fsamp_tr = g_f(fsamp + newLF{fnum}.mu + sparseParam);
        fsamp_tr(find(fsamp_tr<=0.00001)) = 0.00001; %#ok<FNDSB>
        if max(fsamp_tr)>1*max(forceC(fnum,:))
            fsamp_tr = (1*max(forceC(fnum,:))/max(fsamp_tr))*fsamp_tr;
        elseif max(fsamp_tr)<0.8*max(forceC(fnum,:))
            fsamp_tr = (0.8*max(forceC(fnum,:))/max(fsamp_tr))*fsamp_tr;
        end
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
        drawnow();
    end

    % Multiply cascade levels to produce envelopes
    latentGen = bsxfun(@times,lowLevelGen,cascadeMod.fsamp_tr');
    figure(4);clf
    subplot(221)
    plot(HEst2)
    xlim_in = xlim; ylim_in = ylim;
    subplot(223)
    plot(latentGen)
    axis([1 Inf ylim_in])
    drawnow();
    
    latentGen = latentGen * 0.9; % multiply by 0.9 because tNMF always seems to clip
    
    envSynth = WEst2' * latentGen';
    envSynth = makeThisThatLong(envSynth,sig_dur*lfmData.demod_struct.fs)';

    carrier_synth = lfmData.demod_struct.carrier_synth;
    carriers = zeros(size(envSynth));
    for i=1:size(envSynth,1)
        carriers(i,:) = carrier_synth(mod(i-1,size(carrier_synth,1))+1,:);
    end

    sigSynth = sum(carriers'.*envSynth');

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
        audiowrite(strcat('../audio/synth/',fileName,'_TNMF',suffix{suf},'.wav'),sigSynth,lfmData.demod_struct.fs);
    end
    
end

end