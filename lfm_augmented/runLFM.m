function results = runLFM(envInput,entire_sig,guess,do_optim,fixDim,optDim)

% the main script which runs the LFM - input is the amplitude envelopes of
% an audio signal

% Will Wilkinson - Jan 2018

if nargin < 2
    entire_sig = 0;
end
if nargin < 3
    guess = [];
end
if nargin < 4
    do_optim = 1;
end
if nargin < 5
    fixDim = [];
end
if nargin < 6
    optDim = 1:size(envInput.envelopes,1);
    optDim(fixDim) = [];
end

results = struct; % structure for output
results.seed = rng(); % save the seed for reproducability

lfmOpt = envInput.lfmOpt;

%% unpack the settings
analysisRange = lfmOpt.analysisRange; % range of samples to analyse
DS = lfmOpt.DS; % downsample rate

if entire_sig == 1
    inSig = envInput.envelopes([fixDim optDim],1:DS:end);
else
    inSig = envInput.envelopes([fixDim optDim],analysisRange(1):DS:analysisRange(end));
end
if isempty(fixDim) > 0
    optDim = 1:length(optDim);
else
    optDim = fixDim(end)+1:fixDim(end)+length(optDim);
end

% edit the analysis range - exclude periods of silence:
maxT = max(inSig);
silence = find(maxT<lfmOpt.silenceThreshold*max(maxT));
silence_diff = diff(silence);
silence_end = find([silence_diff inf]>1);
silence_start = [1 silence_end(1:end-1)+1];
silence_length = diff([0 silence_end]);
sind = find(silence_length>10);
silence_ind = [silence_start(sind); silence_end(sind)];
silence_ind2 = [];
for i=1:length(sind)
    if silence(silence_ind(2,i)) == size(inSig,2) % if silence ends at the end of the signal then don't adjust
        silence_ind2 = [silence_ind2 silence_ind(1,i)+4:silence_ind(2,i)]; %#ok<AGROW>
    elseif silence(silence_ind(1,i)) == 1 % if silence starts at start of signal then don't adjust
        silence_ind2 = [silence_ind2 silence_ind(1,i):silence_ind(2,i)-4]; %#ok<AGROW>
    else
        silence_ind2 = [silence_ind2 silence_ind(1,i)+4:silence_ind(2,i)-4]; %#ok<AGROW>
    end
end
silence = silence(silence_ind2);
is_sound = 1:length(inSig);
is_sound(silence)=[];
%inSig = inSig_(:,is_sound);

silence_ = find(maxT<lfmOpt.silThrCovReset*max(maxT));
silence_diff_ = diff(silence_);
silence_end_ = find([silence_diff_ inf]>1);
silence_start_ = [1 silence_end_(1:end-1)+1];
silence_length_ = diff([0 silence_end_]);
sind_ = find(silence_length_>10);
silence_ind_ = [silence_start_(sind_); silence_end_(sind_)];
silence_ind2_ = [];
for i=1:length(sind_)
    if silence_(silence_ind_(2,i)) == size(inSig,2) % if silence ends at the end of the signal then don't adjust
        silence_ind2_ = [silence_ind2_ silence_ind_(1,i)+4:silence_ind_(2,i)]; %#ok<AGROW>
    elseif silence_(silence_ind_(1,i)) == 1 % if silence starts at start of signal then don't adjust
        silence_ind2_ = [silence_ind2_ silence_ind_(1,i):silence_ind_(2,i)-4]; %#ok<AGROW>
    else
        silence_ind2_ = [silence_ind2_ silence_ind_(1,i)+4:silence_ind_(2,i)-4]; %#ok<AGROW>
    end
end
silence_ = silence_(silence_ind2_);
is_sound_ = 1:length(inSig);
is_sound_(silence_)=[];

% uncomment if we want to scale the magnitudes
%adjustFactor = 0.1 ./ max(inSig,[],2);
%inSig = bsxfun(@times,inSig,adjustFactor);
%results.adjustFactor = adjustFactor;

nlf = lfmOpt.nlf; % # latent forces
sfactor = lfmOpt.sfactor; % lfm downsample rate (not implemented)
%priorityInference = lfmSettings.priorityInference; % give priority to high amp. outputs (not implemented)
dt = lfmOpt.dt; % time step size
gamma = lfmOpt.gamma; % linearity measure
grng = lfmOpt.grng; % allowable range for gamma
sig_noise = lfmOpt.sig_noise; % assumed measurement noise variance
Dinit = lfmOpt.Dinit; % initial guess for damping
Drng = lfmOpt.Drng; % allowable range of damping values
linit = lfmOpt.linit; % initial guess for lengthscale
lrng = lfmOpt.lrng; % allowable range of lengthscales
Srng = lfmOpt.Srng; % allowable range of sensitivities
resetOpt = lfmOpt.resetOpt; % perform covariance resetting in silent periods during optimisation?
resetFinal = lfmOpt.resetFinal; % perform covariance resetting in silent periods during final filtering and smoothing?
if isfield(lfmOpt,'preset')
    preset = lfmOpt.preset; % period reset time
    preset_it = lfmOpt.preset_it; % number of iterations to perform periodic reset for
else
    preset = 0;
    preset_it = 0;
end
%restrictFunEvals = lfmOpt.restrictFunEvals; % 1 if we restrict number of function evaluations
    MaxFunEvals = lfmOpt.MaxFunEvals;
%restrictIter = lfmOpt.restrictIter; % 1 if we restrict number of iterations
    MaxIter = lfmOpt.MaxIter;
    TolFun = lfmOpt.TolFun;
    
fb_terms = lfmOpt.fb_terms;
fwd_terms = lfmOpt.fwd_terms;
sparsity_fb = lfmOpt.sparsity_fb;
sparsity_fwd = lfmOpt.sparsity_fwd;

if length(sparsity_fb) > fb_terms
    sparsity_fb = sparsity_fb(1:fb_terms);
elseif length(sparsity_fb) < fb_terms
    fb_terms = length(sparsity_fb);
end
sparsity_fb(find(sparsity_fb~=0)) = 1; %#ok<*FNDSB> % must be 0 or 1
if sum(sparsity_fb) == 0
    fb_terms = 0;
    sparsity_fb = [];
end

if length(sparsity_fwd) > fwd_terms
    sparsity_fwd = sparsity_fwd(1:fwd_terms);
elseif length(sparsity_fwd) < fwd_terms
    fwd_terms = length(sparsity_fwd);
end
sparsity_fwd(find(sparsity_fwd~=0)) = 1; % must be 0 or 1
if sum(sparsity_fwd) == 0
    fwd_terms = 0;
    sparsity_fwd = [];
end

% regularisers
lambda_damping = lfmOpt.lambda_damping;
lambda_sensitivity = lfmOpt.lambda_sensitivity;
lambda_feedback = lfmOpt.lambda_feedback;
lambda_lag = lfmOpt.lambda_lag;
%threshold_up_damping = lfmOpt.threshold_up_damping;
%threshold_lo_damping = lfmOpt.threshold_lo_damping;
%threshold_lengthscale = lfmOpt.threshold_lengthscale;

% inital guess for the mean of the latent force at the first time instance:
lf_prior = max(-5.5,log(exp(max(max(inSig(:,2)-inSig(:,1)),0))-1)+1); % if the outputs increase immediately (at time instance 2) then we need input to start immediately


%% Run the LFM
steps = size(inSig,2);

t_min = 1;
t_max = steps;
tgrid = linspace(t_min,t_max,steps);

% Number of outputs
d = size(inSig,1);

% Measurement noise
n_sigma2 = sig_noise^2*ones(d,1);

R = diag(n_sigma2);

% Measurement indices
%msteps = 1:sfactor:steps;
msteps = is_sound; % don't apply update step to silent periods
mind = zeros(1,steps);
mind(msteps) = 1;

msteps_ = is_sound_; % regions for cov. resets
mind_reset = zeros(1,steps);
mind_reset(msteps_) = 1;

% Time between measurements
dty = dt*sfactor;
mysteps = 1:sfactor:steps;
mindy = ones(1,length(mysteps));

% non-linearity
g_param = .1;
% "softplus" rectification function:
g_func = @(f,g_param) log(1+exp(f));
g_df   = @(f,g_param) exp(f)./(1+exp(f));


%% Inference with Gaussian filtering and smoothing:

% Prior variance of the outputs
x_sigma2 = .0001^2;

% General parameters
param = struct;

% Parameters of f
f_func = @lfm_aug_model;
f_param = struct;
f_param.g_func  = g_func;
f_param.g_param = g_param;
f_param.g_df    = g_df;
f_param.x_sigma2 = x_sigma2;
f_param.d = d;

f_param.lf_prior = lf_prior;
f_param.gamma = gamma;
f_param.grng = grng;

% LFM parameters

% Matern model parameters for latent forces
p = 2;
f_param.p = p;
%lengthScale = max(linit,threshold_lengthscale); % (0.5) lengthscale will be optimised
lengthScale = min(max(linit,lrng(1)),lrng(end)); % (0.5) lengthscale will be optimised
sigma = 1.5; % sigma (magnitude) is fixed
f_param.sigma = sigma;

theta_lf = lengthScale*ones(nlf,1);

gj = gamma*ones(d,1); % guess for gamma

Aj = inSig(:,1); % initial conditions
Bj = zeros(d,1);
Dj = Dinit*ones(d,1); % (1.5)
%Sj_sens = 0.5*ones(d,nlf); % 0.5
Sj_sens = bsxfun(@times,max(inSig,[],2),ones(d,nlf));
% This whole loop is simply trying to make a good initial guess for the
% sensitivity parameters by considering the similarity of the observations
if nlf > 1
    shiftMax = floor(1.5*fwd_terms);
    %shiftMax = 1; %%%%%%%%%%%%%%%%%
    mlfGuess = zeros(d,nlf);
    
    [cidx,~,distances] = lfmClustering(inSig,nlf,shiftMax,1,d,100);
    distances = distances(1:d-1)/max(distances(1:d-1));
    
    c_first = zeros(nlf,1);
    c_comp = zeros(nlf,1);
    for j=1:nlf
        c_ = find(cidx==j);
        c_first(j) = c_(1);
        if c_(1) ~= 1
            c_comp(j) = distances(c_(1)-1);
        end
    end
    
    mlfGuess(c_first(1),:) = linspace(0.8,0.2,nlf);
    for j=2:nlf
        mlfGuess(c_first(j),:) = circshift(linspace(0.8,0.2,nlf),j-1,2);
    end
    %{
    mlfGuess(c_first(1),:) = linspace(0.99,0.01,nlf);
    for j=2:nlf
        mlfGuess(c_first(j),:) = circshift(linspace(0.99,0.01,nlf),j-1,2);
    end
    %}
    for i=1:d
        if ismember(i,c_first) == 0
            c_sim = 1 - abs(c_comp - distances(i-1));
            c_sim = c_sim/sum(c_sim); % similarity to previously stored outputs
            for j=1:nlf
                mlfGuess(i,j) = sum(c_sim.*mlfGuess(c_first,j));
            end
        end
    end
    
    Sj_sens = Sj_sens.*mlfGuess;
end
%if isfield(lfmOpt,'mlfGuess')
%    Sj_sens = Sj_sens.*lfmOpt.mlfGuess;
%end
fbn = sum(sparsity_fb); % number of feedback terms to optimise
fwdn = sum(sparsity_fwd); % number of forward terms to optimise
Sj_fb = 10*lfmOpt.fbinit*bsxfun(@times,min(max(inSig,[],2)),ones(d,fbn));%lfmOpt.fbinit*rand(d,fbn);%0.05*rand(d,fbn);
Sj_fwd = 10/nlf*lfmOpt.fwdinit*bsxfun(@times,max(inSig,[],2),ones(d,fwdn*nlf));%lfmOpt.fwdinit*rand(d,fwdn*nlf);%0.025*rand(d,fwdn*nlf);

f_param.A = Aj;
f_param.B = Bj;
f_param.nlf = nlf;
f_param.fbn = fbn;
f_param.fwdn = fwdn;
f_param.fb_terms = fb_terms;
f_param.fwd_terms = fwd_terms;
f_param.sparsity_fb = sparsity_fb;
f_param.sparsity_fwd = sparsity_fwd;

f_param.lambda_damping = lambda_damping;
f_param.lambda_sensitivity = lambda_sensitivity;
f_param.lambda_feedback = lambda_feedback;
f_param.lambda_lag = lambda_lag;
%f_param.threshold_up_damping = threshold_up_damping;
%f_param.threshold_lo_damping = threshold_lo_damping;
%f_param.threshold_lengthscale = threshold_lengthscale;

f_param.fixDim = fixDim;
f_param.optDim = optDim;

if isempty(fixDim) == 0
    if isfield(envInput,'lfmResults')
        [ls_, gamma_, D_, ~, ~, S_sens_, S_fb_, S_fwd_] = unpackTheta(envInput.lfmResults);
        theta_lf = ls_;
        gj(fixDim) = gamma_(fixDim);
        Dj(fixDim) = D_(fixDim);
        Dj(optDim) = mean(D_(fixDim));
        Sj_sens(fixDim,:) = S_sens_(fixDim,:);
        Sj_fb(fixDim,:) = S_fb_(fixDim,find(sparsity_fb==1));
        Sj_fb(optDim,:) = bsxfun(@times,mean(S_fb_(fixDim,find(sparsity_fb==1))),ones(size(Sj_fb(optDim,:))));
        sp_fwd = [];
        for i=1:nlf
            sp_fwd = [sp_fwd sparsity_fwd]; %#ok<AGROW>
        end
        Sj_fwd(fixDim,:) = S_fwd_(fixDim,find(sp_fwd==1));
        Sj_fwd(optDim,:) = bsxfun(@times,quantile(S_fwd_(fixDim,find(sp_fwd==1)),.33),ones(size(Sj_fwd(optDim,:))));
    else
        error('No previous results exist to fix the chosen dimensions')
    end
end

Sj = [Sj_sens(:); Sj_fb(:); Sj_fwd(:)];
%theta = [Dj(:);Sj(:)];
theta = [gj(:);Dj(:);Sj(:)];

% Initial guess for GP lengthscale parameters
%w0_lf = log(theta_lf);
w0_lf = inv_sigmoid(theta_lf,lrng);

% Initial guess for model parameters
%w0    = [log(Dj(:));log(Sj(:))];
%w0    = [inv_sigmoid(Dj(:),Drng);log(Sj(:))];
%w0    = [inv_sigmoid(Dj(:),Drng);inv_sigmoid(Sj(:),Srng)];
w0    = [inv_sigmoid(gj(:),grng);inv_sigmoid(Dj(:),Drng);inv_sigmoid(Sj(:),Srng)];

w0 = [w0_lf;w0];
theta = [theta_lf;theta];

% Parameter prior
theta_prior = cell(1,length(w0));
for i = 1:length(theta)
    theta_prior{i} = prior_t;
    %theta_prior{i}.s2 = .3;
end

lsig_ind = 1:nlf; % lengthscales transformed to defined range using sigmoid
gsig_ind = nlf+1:nlf+length(gj); % gamma transformed to defined range using sigmoid
%Dsig_ind = nlf+1:nlf+length(Dj); % damping transformed to defined range using sigmoid
Dsig_ind = gsig_ind(end)+1:gsig_ind(end)+length(Dj); % damping transformed to defined range using sigmoid
%ltr_ind = Dsig_ind(end)+1:length(theta); % all other parameters must be non-negative so transform them too
Ssig_ind = Dsig_ind(end)+1:length(theta); % sensitivities, fb & fwd transformed to defined range using sigmoid
Sens_sig_ind = Dsig_ind(end)+1:Dsig_ind(end)+length(Sj_sens(:)); % sensitivities
Sfb_ind = Sens_sig_ind(end)+1:Sens_sig_ind(end)+length(Sj_fb(:)); % feedback
Sfwd_ind = Sfb_ind(end)+1:length(theta); % forward

% which parameters to optimise
opt_ind = ones(size(theta));
if isempty(fixDim) == 0
    opt_ind(lsig_ind) = 0;
    opt_ind(gsig_ind(fixDim)) = 0;
    opt_ind(Dsig_ind(fixDim)) = 0;
    for i=0:nlf-1
        opt_ind(Sens_sig_ind(fixDim+i*d)) = 0;
    end
    for i=0:fbn-1
        opt_ind(Sfb_ind(fixDim+i*d)) = 0;
    end
    for i=0:nlf*fwdn-1
        opt_ind(Sfwd_ind(fixDim+i*d)) = 0;
    end
end

if isempty(guess) == 0 && isempty(fixDim) == 1
    if size(guess) == size(w0)
        theta = guess;
        w0 = guess;
        w0(lsig_ind) = inv_sigmoid(w0(lsig_ind),lrng);
        w0(gsig_ind) = inv_sigmoid(w0(gsig_ind),grng);
        w0(Dsig_ind) = inv_sigmoid(w0(Dsig_ind),Drng);
        %w0(ltr_ind) = log(w0(ltr_ind));
        w0(Ssig_ind) = inv_sigmoid(w0(Ssig_ind),Srng);
    else
        error('guess is not the right size')
    end
end

%f_param.opt_ind = opt_ind;
theta_full = theta;
w0_full = w0;
w0 = w0(logical(opt_ind));

% Parameters of q
q_func = @lfm_aug_q;
q_param = struct;
q_param.d = d;
q_param.fb_terms = fb_terms;
q_param.fwd_terms = fwd_terms;

% Setup structures with theta
[param,f_param,q_param] = lfm_aug_set_params(param,f_param,q_param,theta);

M0 = param.M0;
P0 = param.P0;

% Measurement model
h_func = [eye(d) zeros(d,size(M0,1)-d)];
h_param = [];

% store / downsample the data
Y = zeros(d,length(msteps));
obs_ind = 0;
for k=1:steps    
    %if ismember(k,msteps)
        obs_ind = obs_ind + 1;
        Y(:,obs_ind) = inSig(:,k);
    %end
end


%% Inference with fixed parameters
e_param            = struct;
e_param.M0         = M0;
e_param.P0         = P0;
e_param.Y          = Y;
e_param.mind       = mind; %mindy;
e_param.mind_reset = mind_reset; %mindy;
e_param.dt         = dty;
e_param.f_func     = f_func;
e_param.f_param    = f_param;
e_param.q_func     = q_func;
e_param.q_param    = q_param;
e_param.param_func = @lfm_aug_set_params;
e_param.h_func     = h_func;
e_param.h_param    = h_param;
e_param.R          = R;
e_param.t_min      = t_min;
e_param.int_method = @int_cubature;
e_param.lrng       = lrng;
e_param.grng       = grng;
e_param.Drng       = Drng;
e_param.Srng       = Srng;
e_param.lsig_ind   = lsig_ind;
e_param.gsig_ind   = gsig_ind;
e_param.Dsig_ind   = Dsig_ind;
%e_param.ltr_ind    = ltr_ind;
e_param.Ssig_ind   = Ssig_ind;
e_param.resetOpt   = resetOpt;
e_param.resetFinal = resetFinal;
e_param.preset = preset;
e_param.preset_it = preset_it;
e_param.start_ind  = 1;
e_param.ode_func   = @ode45;
e_param.isteps     = 1;
e_param.e_prior = theta_prior;

e_param.d = d;
e_param.nlf = nlf;
e_param.p = p;

e_param.fixDim = fixDim;
e_param.optDim = optDim;
e_param.opt_ind = opt_ind;
e_param.theta_full = theta_full;
e_param.w0_full = w0_full;

e_param.fileName = envInput.demod_struct.fileName;

e_func_gf = @(w) E_sqrt_gf_count(w,e_param);

% Transformed true parameters
ltr_theta = theta;
ltr_theta(lsig_ind) = inv_sigmoid(ltr_theta(lsig_ind),lrng);
ltr_theta(gsig_ind) = inv_sigmoid(ltr_theta(gsig_ind),grng);
ltr_theta(Dsig_ind) = inv_sigmoid(ltr_theta(Dsig_ind),Drng);
%ltr_theta(ltr_ind) = log(ltr_theta(ltr_ind));
ltr_theta(Ssig_ind) = inv_sigmoid(ltr_theta(Ssig_ind),Srng);

%% MAP Optimization of parameter (without gradients)
% Note: the parameter posterior is usually highly non-Gaussian and
% multimodal, so optimization might get stuck to a local mode and even
% fail numerically, depending on simulation and model settings.

tic
fprintf('start time: %s \n',datestr(now))

if do_optim
    opt=optimset('GradObj','off');
    opt=optimset(opt,'TolX', 1e-3);
    opt=optimset(opt,'LargeScale', 'off');
    opt=optimset(opt,'Display', 'iter');
    opt=optimset(opt,'OutputFcn',@outfun);
    %%%%%
    %if restrictFunEvals == 1
        opt.MaxFunEvals = MaxFunEvals;
    %end
    %if restrictIter == 1
        opt.MaxIter = MaxIter;
    %end
        opt.TolFun = TolFun;
    %%%%%
    
    global iterResult; %#ok<*TLEV>
    iterResult = struct;
    iterResult.inSig = inSig;
    iterResult.nlf = nlf;
    iterResult.sparsity_fb = sparsity_fb;
    iterResult.sparsity_fwd = sparsity_fwd;
    iterResult.At = f_param.A;
    iterResult.Bt = f_param.B;
    iterResult.dt = dt;
    %iterResult.gamma = gamma;
    
    global count;
    count = 0;
    global reverseStr;
    reverseStr = '';
    global iterCount;
    iterCount = -1;
    global iterPrev;
    iterPrev = -1;
    warning('off','all'); % too many warnings really slows it down
    w_opt = fminunc(e_func_gf, w0, opt);
    warning('on','all');
    
    w_opt_subset = w_opt;
    w0_full(logical(opt_ind)) = w_opt;
else
    w_opt = w0_full(logical(opt_ind));
    w_opt_subset = w_opt;
end
w_opt = w0_full;

%% filter and smooth with optimised parameters
theta_opt = w_opt;
theta_opt(lsig_ind) = sigmoid(theta_opt(lsig_ind),lrng);
theta_opt(gsig_ind) = sigmoid(theta_opt(gsig_ind),grng);
theta_opt(Dsig_ind) = sigmoid(theta_opt(Dsig_ind),Drng);
%theta_opt(ltr_ind) = exp(theta_opt(ltr_ind));
theta_opt(Ssig_ind) = sigmoid(theta_opt(Ssig_ind),Srng);

theta_opt_subset = theta_opt(logical(opt_ind));

% Interpolation to a finer grid for prettier visualization
isteps = 1;
e_param2 = e_param;
e_param2.dt = dt;
e_param2.mind = mindy;
%e_param2.minds = mind;
e_param2.isteps = isteps;
tgrid2 = linspace(t_min,t_max,steps*isteps);

es_func_gf = @(w) ES_sqrt_gf(w,e_param2);

warning('off','all');

% Estimates with optimized parameters
if do_optim
    [~,MM2,PP2,MS2,PS2] = feval(es_func_gf,w_opt);
    MM1 = MM2;
    PP1 = PP2;
    MS1 = MS2;
    PS1 = PS2;
else
    % Estimates with guessed parameters
    [~,MM1,PP1,MS1,PS1] = feval(es_func_gf,ltr_theta);
    MM2 = MM1;
    PP2 = PP1;
    MS2 = MS1;
    PS2 = PS1;
end
warning('on','all');

toc

%% Plotting estimates of output signals
color1 = [0 0 1];
color2 = [1 0 0];
xx = tgrid'+dt;
xx2 = tgrid2'+dt./isteps;

nisteps = steps*isteps;
H = h_func;
% the predicted output means:
filtOutGuess = H*MM1;
smoothOutGuess = H*MS1;
filtOutOpt = H*MM2;
smoothOutOpt = H*MS2;
% the predicted output variance:
filtOutGuessVar = zeros(d,nisteps);
smoothOutGuessVar = zeros(d,nisteps);
filtOutOptVar = zeros(d,nisteps);
smoothOutOptVar = zeros(d,nisteps);

for i1 = 1:d
    for i = 1:size(PP1,3)
        filtOutGuessVar(i1,i) = H(i1,:)*PP1(:,:,i)*H(i1,:)';
        smoothOutGuessVar(i1,i) = H(i1,:)*PS1(:,:,i)*H(i1,:)';
        filtOutOptVar(i1,i) = H(i1,:)*PP2(:,:,i)*H(i1,:)';
        smoothOutOptVar(i1,i) = H(i1,:)*PS2(:,:,i)*H(i1,:)';
    end
end
figure(3); clf;
for i = 1:d
    subplot(d,2,2*i-1);
    fill([xx2' fliplr(xx2')], [(filtOutOpt(i,:)+1.96*sqrt(abs(filtOutOptVar(i,:)))) ...
        fliplr((filtOutOpt(i,:)-1.96*sqrt(abs(filtOutOptVar(i,:)))))], color1, 'edgecolor',color1);
    hold on
    fill([xx2' fliplr(xx2')], [(filtOutGuess(i,:)+1.96*sqrt(abs(filtOutGuessVar(i,:)))) ...
        fliplr((filtOutGuess(i,:)-1.96*sqrt(abs(filtOutGuessVar(i,:)))))], color2, 'edgecolor',color2); hold on;
    alpha(0.2);
    
    plot(xx2,filtOutOpt(i,:),'color',color1,'LineWidth',2);
    plot(xx2,filtOutGuess(i,:),'color',color2,'LineWidth',2);
    plot(xx,inSig(i,:),'k-','LineWidth',1);
    title(sprintf('Filtered estimate of output %d',i))
    hold off;
    xlim([t_min t_max])
    yl = ylim;
    
    subplot(d,2,2*i);
    fill([xx2' fliplr(xx2')], [(smoothOutOpt(i,:)+1.96*sqrt(abs(smoothOutOptVar(i,:)))) ...
        fliplr((smoothOutOpt(i,:)-1.96*sqrt(abs(smoothOutOptVar(i,:)))))], color1, 'edgecolor',color1);
    hold on
    fill([xx2' fliplr(xx2')], [(smoothOutGuess(i,:)+1.96*sqrt(abs(smoothOutGuessVar(i,:)))) ...
        fliplr((smoothOutGuess(i,:)-1.96*sqrt(abs(smoothOutGuessVar(i,:)))))], color2, 'edgecolor',color2); hold on;
    alpha(0.2);
    
    h1=plot(xx2,smoothOutOpt(i,:),'color',color1,'LineWidth',2);
    h2=plot(xx2,smoothOutGuess(i,:),'color',color2,'LineWidth',2);
    plot(xx,inSig(i,:),'k-','LineWidth',1);
    title(sprintf('Smoothed estimate of output %d',i))
    if i == d
        legend([h1;h2],'Optimised parameters','Guessed parameters')
    end
    hold off;
    xlim([t_min t_max])
    ylim(yl);
    
end
    
%% Plotting of latent forces

Hlf = [zeros(nlf,d) f_param.H];% zeros(1,d*fb_terms+fwd_terms)];
% the predicted latent function means:
filtLFGuess = Hlf*MM1;
smoothLFGuess = Hlf*MS1;
filtLFOpt = Hlf*MM2;
smoothLFOpt = Hlf*MS2;
% the predicted latent function variance:
filtLFGuessVar = zeros(nlf,steps*isteps);
smoothLFGuessVar = zeros(nlf,steps*isteps);
filtLFOptVar = zeros(nlf,steps*isteps);
smoothLFOptVar = zeros(nlf,steps*isteps);

for v = 1:nlf
    for i = 1:size(PP1,3)
        filtLFGuessVar(v,i) = Hlf(v,:)*PP1(:,:,i)*Hlf(v,:)';
        smoothLFGuessVar(v,i) = Hlf(v,:)*PS1(:,:,i)*Hlf(v,:)';
        filtLFOptVar(v,i) = Hlf(v,:)*PP2(:,:,i)*Hlf(v,:)';
        smoothLFOptVar(v,i) = Hlf(v,:)*PS2(:,:,i)*Hlf(v,:)';
    end
end

figure(4); clf;
subplot(2,1,1);
fill([xx2' fliplr(xx2')], [(filtLFOpt+1.96*sqrt(abs(filtLFOptVar))) ...
    fliplr((filtLFOpt-1.96*sqrt(abs(filtLFOptVar))))], color1, 'edgecolor',color1);
hold on
fill([xx2' fliplr(xx2')], [(filtLFGuess+1.96*sqrt(abs(filtLFGuessVar))) ...
    fliplr((filtLFGuess-1.96*sqrt(abs(filtLFGuessVar))))], color2, 'edgecolor',color2); hold on;
alpha(0.2);

plot(xx2,filtLFOpt,'color',color1,'LineWidth',2);
plot(xx2,filtLFGuess,'color',color2,'LineWidth',2);
title('Filtered estimate of force')
hold off;
xlim([t_min t_max])
yl = ylim;

subplot(2,1,2);
fill([xx2' fliplr(xx2')], [(smoothLFOpt+1.96*sqrt(abs(smoothLFOptVar))) ...
    fliplr((smoothLFOpt-1.96*sqrt(abs(smoothLFOptVar))))], color1, 'edgecolor',color1);
hold on
fill([xx2' fliplr(xx2')], [(smoothLFGuess+1.96*sqrt(abs(smoothLFGuessVar))) ...
    fliplr((smoothLFGuess-1.96*sqrt(abs(smoothLFGuessVar))))], color2, 'edgecolor',color2); hold on;
alpha(0.2);

h1=plot(xx2,smoothLFOpt,'color',color1,'LineWidth',2);
h2=plot(xx2,smoothLFGuess,'color',color2,'LineWidth',2);
title('Smoothed estimate of force')
legend([h1(1);h2(1)],'Optimised parameters','Guessed parameters')
hold off;
xlim([t_min t_max])
ylim(yl);


%% Save results
results.theta_opt = theta_opt;
results.Y = Y;

results.inSig = inSig;

results.OutSMAP = smoothOutOpt;
results.VarS = smoothOutOptVar;
results.OutFMAP = filtOutOpt;
results.VarF = filtOutOptVar;

results.ForceSMAP = smoothLFOpt;
results.ForceFMAP = filtLFOpt;
results.uVarS = smoothLFOptVar;
results.uVarF = filtLFOptVar;
results.gS = log(1+exp(smoothLFOpt));
results.gF = log(1+exp(filtLFOpt));

results.nlf = nlf;
results.p = p;
results.sigma = sigma;

results.dt = dt;
results.dty = dty;
results.sfactor = sfactor;
%results.gamma = gamma;

results.At = Aj;
results.Bt = Bj;
results.lsig_ind = lsig_ind;
results.gsig_ind = gsig_ind;
results.Dsig_ind = Dsig_ind;
%results.ltr_ind = ltr_ind;
results.Ssig_ind = Ssig_ind;

results.fb_terms = fb_terms;
results.fwd_terms = fwd_terms;
results.sparsity_fb = sparsity_fb;
results.sparsity_fwd = sparsity_fwd;

results.e_param = e_param;
results.e_param2 = e_param2;

results.w_opt_subset = w_opt_subset;
results.theta_opt_subset = theta_opt_subset;

end