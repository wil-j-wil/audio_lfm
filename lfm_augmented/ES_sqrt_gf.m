% ES_SQRT_GF Log likelihood and state estimates with sqrt Gaussian filter
% 
% Syntax:
%   [E,MM,PP,MS,PS,MPR,PPR,GG] = ES_sqrt_gf(theta,param)
%
% In:
%   theta - Parameters as dx1 vector
%   param - Parameter structure (see below for details)
%      
% Out:
%   E - Negative log likelihood
%   MM - Filtered state means
%   PP - Filtered state covariances
%   MS - Smoothed state means
%   PS - Smoothed state covariances
%   MPR - Predictive state means
%   PPR - Predictive state covariance
%   GG  - Smoother gains
% 
% Description:
%  Calculates the negative log likelihood of data and state estimates with
%  a partially observed SDE model 
% 
%  dx(t)/dt = f(x(t),t)   + L(x(t),t) w(t)
%       y_k = h_k(x(t_k)) + r_k, r_k ~ N(0,R_k), k = 1,...,T
% 
%  using a square-root form continuous-discrete Gaussian filter. 
%
% Initially based on code by Jouni Hartikainen, Simo Särkkä
%
% Modified for audio LFM work by Will Wilkinson - Jan 2018

function [E,MM,PP,MS,PS,MPR,PPR,GG] = ES_sqrt_gf(theta, param)

% Extract parameters
Y          = param.Y;           % Measurements
mind       = param.mind;        % Measurement index vector
mind_reset = param.mind_reset;  % Regions where cov. should be reset
reset_cov  = param.resetFinal;  % Flag whether to perform covariance resetting
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
isteps     = param.isteps;      % Number of Runge-Kutta steps

d          = param.d;           % Number of outputs
nlf        = param.nlf;         % Number of latent forces
p1         = param.p + 1;       % GP order (+1), i.e. size of SDE state space
fb_terms   = f_param.fb_terms;  % Number of feedback terms
fwd_terms  = f_param.fwd_terms; % Number of forward terms
lf_prior   = f_param.lf_prior;  % inital guess for the mean of the latent force at the first time instan

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
end

% Parameters of ODE solver
if isfield(param,'ode_opt')
    opt = param.ode_opt;
else
    opt = odeset;
end

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

dtrk = dt/isteps;

% Transform the parameters
theta(lsig_ind) = sigmoid(theta(lsig_ind),lrng);
theta(gsig_ind) = sigmoid(theta(gsig_ind),grng);
actual_g  = theta(gsig_ind);
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
steps = length(mind);

% Form X and W for the chosen sigma point rule
[XI,WM,WC] = feval(int_method,n,int_param);

% Parameters for the moment derivative function to be integrated
ode_param = {XI,WM,WC,f_func,f_param,Q,t_min};

% Space for estimates
AA = zeros(n,n,steps*isteps);
PP = zeros(n,n,steps*isteps);
MM = zeros(n,steps*isteps);
CC = zeros(n,n,steps*isteps);
AAR = zeros(n,n,steps*isteps);
PPR = zeros(n,n,steps*isteps);
MPR = zeros(n,steps*isteps);
EE = zeros(1,steps);

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

% Measurement counter
mc = 1;

n_l = d*fb_terms+nlf*fwd_terms; % number of lagged terms
m_l = zeros(n_l,1); % prior lagged mean
A_l = zeros(n_l); % prior lagged covariance
%global envCount;
% append lagged terms
m = [m; m_l];
A = blkdiag(A,A_l);
A_prior = A; % prior covariance
rcount = 0;
k_reset = [];
%reset_flag = 0;
% Filtering
for k=1:steps
    %k
    if reset_cov == 1
        %k____ = k
        if mind_reset(k) == 0
            rcount = rcount + 1;
        else
            rcount = 0;
        end
        %{
        if (k > 1 && mind_reset(k) ~= mind_reset(k-1)) || rcount > 75
            A = 0.01*A_prior + 0.99*A; % reset the cholesky & covariance to prior at start/end of silence
            k_reset = [k_reset; k]; %#ok<AGROW>
            rcount = 0;
        end
        %}
        if (k > 1 && mind_reset(k) ~= mind_reset(k-1))
            if mind_reset(k) == 0
                opt.RelTol = 0.1;
                opt
            else
                opt.RelTol = [];
                opt
            end
        end
        if k==1 && mind_reset(k) == 0
            opt.RelTol = 0.1;
            opt
        end
        
    end    
    
    m_ = m;
    if k > start_ind        
        t = t_min+dt*(k-1);
        tt = t;
        for i = 1:isteps
            mAC = [m A A*A'];
            %if ismember(k,k_reset)
            %    mAC = [m A C_];
            %end
            vmAC = mAC2vec(mAC,n+n_l);
            %m'
            %pause
            Tspan = [tt tt+dtrk];
            tt = tt + dtrk;
            % Solve the moment ODEs
            [sol,y] = feval(ode_func,@(t,y) myode(t,y,@lfm_aug_dmAC_sqrt,ode_param),Tspan,vmAC,opt); %#ok<ASGLU>
            vmAC = y(end,:)';
            [mAC,m,A,C] = vec2mAC(vmAC,n+n_l); %#ok<ASGLU>
            % Save the moments
            MPR(:,i+(k-1)*isteps)   = m(1:n);
            AAR(:,:,i+(k-1)*isteps) = A(1:n,1:n);
            PPR(:,:,i+(k-1)*isteps) = A(1:n,1:n)*A(1:n,1:n)';
            CC(:,:,i+(k-1)*isteps)  = C(1:n,1:n);
            MM(:,i+(k-1)*isteps)    = m(1:n);
            AA(:,:,i+(k-1)*isteps)  = A(1:n,1:n);
            PP(:,:,i+(k-1)*isteps)  = A(1:n,1:n)*A(1:n,1:n)';
        end
    end
    
    if mind(k) == 1
        % Linear update
        [m,P,A,W,LLH] = srkf_update(m(1:n),A(1:n,1:n),Y(:,mc),h_func,chol(R)'); %#ok<ASGLU>
        A = blkdiag(A,A_l);
        mc = mc + 1;
        EE(k) = LLH;
    end

    MM(:,k*isteps)   = m(1:n);
    AA(:,:,k*isteps) = A(1:n,1:n);
    PP(:,:,k*isteps) = A(1:n,1:n)*A(1:n,1:n)';
    
    % correct the non-Markovian terms (we need them to be exact)
    if fb_terms > 0
        m(fb1) = m_(1:d); % first fb terms = previous predicted outputs
        m(fb2) = m_(fb2-1); % bucket brigade fb terms
    end
    if fwd_terms > 0
        m(fwd1) = log(1+exp(m_(gpind))); % first fwd terms = previous GP terms
        m(fwd2) = m_(fwd2-1); % bucket brigade fwd terms
    end
    
    max(max(isnan(m)))
end

E = -sum(EE) + E_prior;


%% Smoothing
MS = MM;
PS = zeros(size(AA));
PS(:,:,end) = AA(:,:,end)*AA(:,:,end)';
ms = MM(:,end);
Ps = AA(:,:,end)*AA(:,:,end)';
GG = zeros(n,n,isteps*steps);
for k=isteps*steps-1:-1:1
    if (mind_reset(k) == 1 && ismember(k+1,k_reset) == 0) || reset_cov == 0
    %if (mind_reset(k) == 1 && ismember(k,k_reset) == 0) || reset_cov == 0

        MPr = MPR(:,k+1);
        PPr = AAR(:,:,k+1)*AAR(:,:,k+1)';
        G = (CC(:,:,k+1)/AAR(:,:,k+1)')/AAR(:,:,k+1);
        m = MM(:,k);
        P = AA(:,:,k)*AA(:,:,k)';
        ms = m + G * (ms - MPr);
        Ps = P + G * (Ps - PPr) * G';
        MS(:,k) = ms;
        PS(:,:,k) = Ps;
        GG(:,:,k) = G;
    else % during silent periods, just use the filtering solution
        PS(:,:,k) = AA(:,:,k)*AA(:,:,k)';
        ms = MM(:,k);
        Ps = AA(:,:,k)*AA(:,:,k)';
    end
end


