% Edited version of TF_SET_PARAMS Parameter setup function of transcription
% factor model from LFM toolbox
% 
% Syntax:
%   [param f_param q_param] = tf_set_params(param,f_param,q_param,theta)
%
% In:
%   param   - General parameter structure
%   f_param - Parameters of dynamic model
%   q_param - Parameters of dispersion model
%   theta   - Vector of parameters (alpha,beta,sigma)
% Out:
%   param   - Updated general parameter structure
%   f_param - Updated dynamic model parameter structure
%   q_param - Updated dispersion model parameter structure
%
% Description:
%   Parameter setup function for transcriptional factor model, described by
%   the equations
%   
%   dx_i(t)/dt = B_i + \sum_{r=1}^R S_{i,r} g_i(u_r(t)) - D_i x(t),
% 
%   for outputs x_i(t), i=1,...,N and latent forces u_r(t), r=1,...,R,
%   which modelled with LTI SDE models specified by the user in param
%   structure. The non-linear function g_i(.) is also specified in param.
%
% Initially based on code by Jouni Hartikainen, Simo Särkkä
%
% Modified for audio LFM work by Will Wilkinson - Jan 2018

function [param, f_param, q_param] = lfm_aug_set_params(param,f_param,q_param,theta)
% Latent force parameters

% Matern model:
p = f_param.p;
nlf = f_param.nlf;

ls = theta(1:nlf);
%sigma = theta(2);
sigma = f_param.sigma;
%sigma_ = f_param.sigma;

%F = model.F;
%H = model.H;
%q = model.Qc;
%M0gp = model.M0;
%P0gp = model.P0;
F = [];
H = [];
q = [];
M0gp = [];
P0gp = [];
for v=1:nlf
    model = matern_model([ls(v) sigma],p);
    F = blkdiag(F,model.F);
    H = blkdiag(H,model.H);
    q = blkdiag(q,model.Qc);
    M0gp = [M0gp; model.M0]; %#ok<*AGROW>
    P0gp = blkdiag(P0gp,model.P0);
end

q_param.Qgp = q;

f_param.F = F;
f_param.H = H;

% Output parameters
d = f_param.d;

fb_terms = f_param.fb_terms;
fwd_terms = f_param.fwd_terms;
fbn = f_param.fbn;
fwdn = f_param.fwdn;

%ind = 2*nlf+1:2*nlf+d; % first 2 params are the GP hyperparameters (manual for now)
ind = nlf+1:nlf+d; % index to gamma (first nlf params are the GP hyperparameters)
gj = theta(ind); % gamma
ind = ind(end)+1:ind(end)+d; % index to damping
%ind_ = 2:1+d;
%Aj = theta(ind);
%ind = ind+d;
%Bj = theta(ind);
%ind = ind+d;
Dj = theta(ind);
ind = ind(end)+1:ind(end)+d*nlf;
Sj_sens = theta(ind);
if fbn > 0
    ind = ind(end)+1:ind(end)+d*fbn;
    Sj_fb = theta(ind);
else
    Sj_fb = zeros(d,0);
end
ind = ind(end)+1:ind(end)+d*fwdn*nlf;
Sj_fwd = theta(ind);

%f_param.A = Aj;
%f_param.B = Bj;
f_param.gamma = gj;
f_param.D = Dj;
f_param.S_sens = Sj_sens;
f_param.S_fb = Sj_fb;
f_param.S_fwd = Sj_fwd;

x_sigma2 = f_param.x_sigma2;

%param.M0 = [Aj+Bj./Dj;M0gp];
%param.M0 = [Aj;M0gp];

% fb
A_fb = [];
for i=1:fb_terms
    A_fb = [A_fb;f_param.A];
end
% fwd
A_fwd = [];
for i=1:fwd_terms
    A_fwd = [A_fwd;H*M0gp];
end

%param.M0 = [f_param.A;M0gp;A_fb;A_fwd];
%param.P0 = blkdiag(x_sigma2*eye(d),P0gp);
%param.P0 = blkdiag(x_sigma2*eye(d),P0gp,x_sigma2*eye(d*fb_terms+fwd_terms));

param.M0 = [f_param.A;M0gp];
param.M0b = [A_fb;A_fwd];
param.P0 = blkdiag(x_sigma2*eye(d),P0gp);
param.P0b = blkdiag(x_sigma2*eye(d*fb_terms+fwd_terms));

% if covariance is not positive semi-definite then find the nearest matrix
% that is
[~,psdCheck] = chol(param.P0);
if psdCheck > 0
    param.P0 = nearestSPD(param.P0);
end

end