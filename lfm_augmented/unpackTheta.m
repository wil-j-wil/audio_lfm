function [ls, gamma, D, S, S_, S_sens, S_fb_, S_fwd_] = unpackTheta(results)

% unpack the model parameters stored in the theta vector

% Will Wilkinson - Jan 2018

% ls = lengthscale (GP hyperparameter)
% sigma = magnitude (GP hyperparameter)
% gamma = linearity measure
% D = Damping
% S = Sensitivity
% S_ = Sensitivity (reformatted with sparsity)
% S_sens = first order sensitivity
% S_fb_ = feedback
% S_fwd_ = forward / distributed lag / higher-order sensitivity

N = size(results.inSig,1);
nlf = results.nlf;
theta = results.theta_opt;
sparsity_fb = results.sparsity_fb;
sparsity_fwd = results.sparsity_fwd;
fbn = sum(sparsity_fb);
fwdn = sum(sparsity_fwd);
fb_terms = length(sparsity_fb);
fwd_terms = length(sparsity_fwd);

% GP hyperparameters
%ls = theta(1:2:2*nlf);
ls = theta(1:nlf);
%sigma = theta(2:2:2*nlf+1);
gind = nlf+1:nlf+N; % theta contains 1 hyperparameter per GP/LF
%Dind = 2*nlf+1:2*nlf+N; % theta contains 2 hyperparameters per GP/LF
Dind = gind(end)+1:gind(end)+N; % index to damping
Ssind = Dind(end)+1:Dind(end)+N*nlf;
%{
Sfbind = Ssind(end)+1:Ssind(end)+N*fbn;
if fwdn > 0
    Sfwdind = Sfbind(end)+1:Sfbind(end)+N*fwdn*nlf;
else
    Sfwdind = zeros(0,1);
end
%}
if fbn > 0
    Sfbind = Ssind(end)+1:Ssind(end)+N*fbn;
    if fwdn > 0
        Sfwdind = Sfbind(end)+1:Sfbind(end)+N*nlf*fwdn;
    else
        Sfwdind = zeros(0,1);
    end
else
    Sfbind = zeros(0,1);
    if fwdn > 0
        Sfwdind = Ssind(end)+1:Ssind(end)+N*nlf*fwdn;
    else
        Sfwdind = zeros(0,1);
    end
end


gamma = theta(gind);
D = theta(Dind);
S = [theta(Ssind); theta(Sfbind); theta(Sfwdind)];

S_sens = reshape(theta(Ssind),N,nlf);
S_fb_ = zeros(N,fb_terms);
S_fb_(:,find(sparsity_fb==1)) = reshape(theta(Sfbind),N,fbn); %#ok<*FNDSB>
S_fwd_ = zeros(N,fwd_terms*nlf);
if fwdn > 0
    spind = [];%find(sparsity_fwd==1);
    for v=0:nlf-1
        spind = [spind v*fwd_terms+find(sparsity_fwd==1)]; %#ok<*AGROW>
    end
    S_fwd_(:,spind) = reshape(theta(Sfwdind),N,fwdn*nlf);
end

S_ = [S_sens S_fb_ S_fwd_];

end