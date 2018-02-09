% E_SQRT_GF Negative log likelihood with sqrt Gaussian filter
% 
% Syntax:
%  E = E_sqrt_gf(theta,param)
%
% In:
%   theta - Parameters as dx1 vector
%   param - Parameter structure (see below for details)
%      
% Out:
%   E - Negative log likelihood
% 
% Description:
%  Calculates the negative log likelihood of data with a partially observed
%  SDE model 
% 
%  dx(t)/dt = f(x(t),t)   + L(x(t),t) w(t)
%       y_k = h_k(x(t_k)) + r_k, r_k ~ N(0,R_k), k = 1,...,T
% 
%  using a square-root form continuous-discrete Gaussian filter. 
%
% Initially based on code by Jouni Hartikainen, Simo Särkkä
%
% Modified for audio LFM work by Will Wilkinson - Jan 2018

function [E] = E_sqrt_gf_count(theta, param)

Y          = param.Y;           % Measurements
mind       = param.mind;        % Measurement index vector
mind_reset = param.mind_reset;  % Regions where cov. should be reset
reset_cov  = param.resetOpt;    % Flag whether to perform covariance resetting
preset     = param.preset;      % period reset?
preset_it  = param.preset_it;   % period reset iterations
f_func     = param.f_func;      % Dynamic model function
f_param    = param.f_param;     % Parameters of f 
q_func     = param.q_func;      % Dispersion function
q_param    = param.q_param;     % Parameters of q   
param_func = param.param_func;  % Function for setting parameter vectors w. theta
int_method = param.int_method;  % Moment integration method
h_func     = param.h_func;      % Measurement model function
%h_param    = param.h_param;     % Parameters of h
R          = param.R;           % Measurement noise covariance (matrix or function returning matrix)
dt         = param.dt;          % Time step between measurements   
t_min      = param.t_min;       % Time at initial time step
lrng       = param.lrng;        % range of allowable values for lengthscales
grng       = param.grng;        % range of allowable values for gamma
Drng       = param.Drng;        % range of allowable values for damping
Srng       = param.Srng;        % range of allowable values for sensitivity
lsig_ind   = param.lsig_ind;    % Indicator vector for sigmoid-transformed lengthscale parameters
gsig_ind   = param.gsig_ind;    % Indicator vector for sigmoid-transformed gamma (linearity) parameters
Dsig_ind   = param.Dsig_ind;    % Indicator vector for sigmoid-transformed damping parameters
%ltr_ind    = param.ltr_ind;     % Indicator vector for log-transformed parameters
Ssig_ind   = param.Ssig_ind;    % Indicator vector for sigmoid-transformed sensitivity parameters

d          = param.d;           % Number of outputs
nlf        = param.nlf;         % Number of latent forces
p1         = param.p + 1;       % GP order (+1), i.e. size of SDE state space
fb_terms   = f_param.fb_terms;  % Number of feedback terms
fwd_terms  = f_param.fwd_terms; % Number of forward terms
lf_prior   = f_param.lf_prior;  % inital guess for the mean of the latent force at the first time instance

theta_full = param.w0_full;     % full theta vector including fixed parameters
opt_ind    = param.opt_ind;     % index to parameters being optimised
fileName   = param.fileName;

theta_ = theta_full;
theta_(logical(opt_ind)) = theta;
theta = theta_;

% Parameters of R if it's a function
if isfield(param,'n_param')
    n_param = param.n_param;
else
    n_param = [];
end

if ~isnumeric(R)
    R = R(theta,n_param);
end

% Function used to solve ODEs of m and P wrt time
if isfield(param,'ode_func')
    ode_func = param.ode_func;
else
    % This is usually the best choice
    ode_func = @ode45;
    %ode_func = @ode15s; % WJW - Mar 2017
end

% Parameters of ODE solver
if isfield(param,'ode_opt')
    opt = param.ode_opt;
else
    opt = odeset;
end
%opt.NonNegative = 1:6;
%opt.MaxOrder = 5; % Could set to 1 to make faster. Not sure how it would affect accuracy. WJW - Mar 2017

if isfield(param,'start_ind')
    start_ind = param.start_ind;
else
    start_ind = 0;
end

% Parameters of moment integration function
if isfield(param,'int_param')
    int_param = param.int_param;
else
    int_param = {};
end

% Transform the parameters
theta(lsig_ind) = sigmoid(theta(lsig_ind),lrng);
theta(gsig_ind) = sigmoid(theta(gsig_ind),grng);
theta(Dsig_ind) = sigmoid(theta(Dsig_ind),Drng);
%theta(ltr_ind) = exp(theta(ltr_ind));
theta(Ssig_ind) = sigmoid(theta(Ssig_ind),Srng);

% Add prior contribution
E_prior = 0;
if isfield(param,'e_prior')
    for i = 1:length(theta)
        E_prior = E_prior - real(feval(param.e_prior{i}.fh.lp,theta(i),param.e_prior{i}));   
    end
end

% Set theta to param, f_param and q_param
[param, f_param, q_param] = feval(param_func,param,f_param,q_param,theta);

% Dispersion matrix
Q = feval(q_func,[],[],q_param);

% This notation is used with Gaussian filters
Q = Q*Q';

% Prior mean and covariance
M0 = param.M0;
P0 = param.P0;

m = M0;
for v=0:nlf-1
    m(d+1+v*p1) = lf_prior;
end
P = P0;
A = chol(P)';

n = size(m,1);
%steps = length(mind);
steps = size(Y,2);

% Space for conditional measurement likelihoods
EE = zeros(1,steps);

% Form X and W for the chosen sigma point rule
[XI,WM,WC] = feval(int_method,n,int_param);

% Parameters for the moment derivative function to be integrated
ode_param = {XI,WM,WC,f_func,f_param,Q,t_min};

% Regularisation
theta_damp = f_param.D;
theta_sens = f_param.S_sens;
theta_fb = f_param.S_fb;
theta_fwd = f_param.S_fwd;
%theta_l = theta(1:nlf); % lengthscales

%theta_damp_inv = 1./theta_damp; % encourage high damping
%rgls_damp = f_param.lambda_damping*sqrt(theta_damp_inv'*theta_damp_inv); % regularise inverse damping with L2 norm
rgls_damp = f_param.lambda_damping*sqrt(theta_damp'*theta_damp); % regularise inverse damping with L2 norm
rgls_sens = f_param.lambda_sensitivity*sqrt(theta_sens'*theta_sens); % regularise sesnitivity with L2 norm
if isempty(theta_fb)
    rgls_fb = 0;
else
    rgls_fb = f_param.lambda_feedback*sum(abs(theta_fb)); % regularise feedback with L1 norm
    %rgls_fb = f_param.lambda_feedback*sqrt(theta_fb'*theta_fb); % regularise feedback with L2 norm
end
if isempty(theta_fwd)
    rgls_fwd = 0;
else
    rgls_fwd = f_param.lambda_lag*sum(abs(theta_fwd)); % regularise distributed lag with L1 norm
    %rgls_fwd = f_param.lambda_lag*sqrt(theta_fwd'*theta_fwd); % regularise distributed lag with L2 norm
end

% Measurement counter
%mm = [];
MM = [];
mc = 1;

gpind = d+1:p1:d+nlf*p1; % index to GP terms

fbind = d+nlf*p1+1:d+nlf*p1+(d*fb_terms); % index to all feedback terms
if fb_terms > 0
    fb1 = fbind(1):fb_terms:fbind(end); % index to the first feedback terms
    fb2 = fbind; % index to the remaining feedback terms
    for l=fb1
        fb2 = fb2(fb2~=l);
    end
%if fb_terms > 0
    fwdind = fbind(end)+1:fbind(end)+nlf*fwd_terms; % index to all forward terms
else
    fwdind = d+nlf*p1+1:d+nlf*p1+nlf*fwd_terms;
end
if fwd_terms > 0
    fwd1 = fwdind(1):fwd_terms:fwdind(end); % index to the first forward terms
    fwd2 = fwdind; % index to the remaining forward terms
    for l=fwd1
        fwd2 = fwd2(fwd2~=l);
    end
end

n_l = d*fb_terms+nlf*fwd_terms; % number of lagged terms
m_l = zeros(n_l,1); % prior lagged mean
A_l = zeros(n_l); % prior lagged covariance

% append lagged terms
m = [m; m_l];
A = blkdiag(A,A_l);
A_prior = A; % prior covariance
rcount = 0;
pcount = 0;
global iterCount
if iterCount > preset_it
    preset = 0; % don't periodically reset after first iteration
end
%k_reset = [];
%reset_flag = 0;
% Filtering
for k=1:steps
    %if k > 1 && mind_reset(k) ~= mind_reset(k-1)
    %    A = A_prior; % reset the covariance to prior at start/end of silence
    %end
    pcount = pcount + 1;
    if preset > 0 && pcount > preset
        A = 0.01*A_prior + 0.99*A; % reset the cholesky & covariance to prior at start/end of silence
        pcount = 0;
    end
    
    if reset_cov == 1
        if mind_reset(k) == 0
            rcount = rcount + 1;
        else
            rcount = 0;
        end
        if (k > 1 && mind_reset(k) ~= mind_reset(k-1)) || rcount > 75
            A = 0.01*A_prior + 0.99*A; % reset the cholesky & covariance to prior at start/end of silence
            %k_reset = [k_reset; k]; %#ok<AGROW>
            rcount = 0;
        end
        
    end
    
    m_ = m;
    if mind(k) == 1
        if k > start_ind
            % append lagged terms
            %m_temp = [m; m_l];
            %A_temp = blkdiag(A,A_l);
            %mA = [m_temp A_temp];
            mA = [m A];
            vmA = mA2vec(mA,n+n_l);
            t = t_min+dt*(k-1);
            Tspan = [t t+dt];
            [sol,y] = feval(ode_func,@(t,y) myode(t,y,@lfm_aug_dmA_sqrt,ode_param),Tspan,vmA,opt); %#ok<ASGLU>
            vmA = y(end,:)';
            [mA,m,A] = vec2mA(vmA,n+n_l); %#ok<ASGLU>
        end
    else
            %for j=0:nlf-1
                %m(gpind+j*p1) = min(log(exp(0.001)-1) , min(min(MM(gpind,:)))); % some value that is close to zero when passed through the nonlinearity
            %end
        
            %m(gpind) = min(log(exp(0.001)-1) , min(min(MM(gpind,:)))); % some value that is close to zero when passed through the nonlinearity
        m(gpind) = log(exp(max(Y(:,k)))-1); % some value that is close to zero when passed through the nonlinearity
    end
    %mm = [mm m];
    %if mind(k) == 1
        % Linear update (only the outputs and GP terms)
        [m,P,A,W,LLH] = srkf_update(m(1:n),A(1:n,1:n),Y(:,mc),h_func,chol(R)'); %#ok<ASGLU>
        A = blkdiag(A,A_l);
        mc = mc + 1;
        if mind(k) == 1
            EE(k) = LLH;
        end
        
        %if k > start_ind
            m = [m; m_(n+1:n+n_l)]; %#ok<AGROW>
            % correct the non-Markovian terms (we need them to be exact)
            if fb_terms > 0
                m(fb1) = m_(1:d); % first fb terms = previous predicted outputs
                m(fb2) = m_(fb2-1); % bucket brigade fb terms
            end
            if fwd_terms > 0
                m(fwd1) = log(1+exp(m_(gpind))); % first fwd terms = previous GP terms
                m(fwd2) = m_(fwd2-1); % bucket brigade fwd terms
            end
        %end
        
    %end
    MM = [MM m]; %#ok<AGROW>
    % UNCOMMENT THIS k TO GET A STEP BY STEP COUNTER TO SEE WHERE THE
    % FILTERING IS SLOW / BREAKS DOWN:
    %k
end

E = -sum(EE) + E_prior + rgls_damp + rgls_sens + rgls_fb + rgls_fwd;% + thr_up_damp + thr_lo_damp + thr_l; % likelihood + regularisation

Hlf = [zeros(nlf,d) f_param.H];
global iterResult
iterResult.ForceFMAP = Hlf*MM(1:size(Hlf,2),:);
iterResult.theta = theta;
iterResult.lrng = lrng;
iterResult.grng = grng;
iterResult.Drng = Drng;
iterResult.Srng = Srng;
iterResult.lsig_ind = lsig_ind;
iterResult.gsig_ind = gsig_ind;
iterResult.Dsig_ind = Dsig_ind;
iterResult.Ssig_ind = Ssig_ind;
%iterResult.ltr_ind = ltr_ind;
iterResult.fileName = fileName;

global count
count = count + 1;
global reverseStr
%global iterCount
global iterPrev
if iterCount == -1
    msg=sprintf('Func-count %d',count);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
elseif iterCount == iterPrev
    msg=sprintf('Func-count %d \n',count);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
else
    reverseStr = '';
    iterPrev = iterCount;
end