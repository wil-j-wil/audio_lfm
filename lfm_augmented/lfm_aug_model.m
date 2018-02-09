% Edited version of TF_F from the LFM toolbox. Now includes
% gamma parameter and extra feedback and lag terms

% TF_F Transcription factor model (1st order non-linear LFM) dynamics
% 
% Syntax:
%   f = tf_f(x,t,param)
%
% In:
%   x     - State vector
%   t     - Time instance
%   param - Model parameters (see below for details)    
% Out:
%   f     - Jacobian of f(x(t),t) wrt x
%
% Description:
%   Dynamics of transcriptional factor, described by the equations
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

function f = lfm_aug_model(x,t,param)
    
    % Extract needed parameters
    A = param.A;
    B = param.B;
    D = param.D;
    S_sens = param.S_sens;
    S_fb = param.S_fb;
    S_fwd = param.S_fwd;
    d = param.d; % number of outputs
    nlf = param.nlf; % number of latent force
    F = param.F;
    H = param.H;
    g = param.g_func;
    g_param = param.g_param;
    
    gamma = param.gamma;
    
    fb_terms = param.fb_terms;
    fwd_terms = param.fwd_terms;
    sparsity_fb = param.sparsity_fb;
    sparsity_fwd = param.sparsity_fwd;

    % Storage for time derivative
    f = zeros(size(x));   
    % Evaluate transformed forces
    gf = feval(g,H*x(d+1:d+size(H,2),:),g_param);
    
    gp_ind = d+1:d+size(H,2);
    if fb_terms > 0
        fb_ind = gp_ind(end)+1:gp_ind(end)+d*fb_terms;
        if fwd_terms > 0
            fwd_ind = fb_ind(end)+1:fb_ind(end)+nlf*fwd_terms;
        else
            fwd_ind = zeros(0,1);
        end
    else
        fb_ind = zeros(0,1);
        if fwd_terms > 0
            fwd_ind = gp_ind(end)+1:gp_ind(end)+nlf*fwd_terms;
        else
            fwd_ind = zeros(0,1);
        end
    end
    
    
    S_sens = reshape(S_sens,d,nlf); % reshape the sensitivity parameters
    
    if fb_terms > 0
        S_fb_ = zeros(fb_terms*d,1);
        ind = 1:d:length(S_fb);
        fb1=find(sparsity_fb);
        for i=0:d-1
            S_fb_(fb1+i*fb_terms) = S_fb(ind);
            ind = ind + 1;
        end
    end
    
    if fwd_terms > 0
        S_fwd_ = [];
        jind = 1:d;
        for i=1:nlf
            for j=1:fwd_terms
                if sparsity_fwd(j)
                    S_fwd_ = [S_fwd_; S_fwd(jind)];
                    jind = jind + d;
                else
                    S_fwd_ = [S_fwd_; zeros(d,1)];
                end
            end
        end
        S_fwd_=reshape(S_fwd_,d,fwd_terms*nlf);
    end
    
    % output dynamics
    if isempty(fwd_ind) && isempty(fb_ind) % no high-order terms
        for i1 = 1:size(x,2)
            f(1:d,i1) = S_sens*gf(:,i1) - D.*real(x(1:d,i1).^(gamma));
        end
    elseif isempty(fwd_ind) % feedback terms only
        for i1 = 1:size(x,2)
            f(1:d,i1) = S_sens*gf(:,i1) + sum(reshape(S_fb_.*x(fb_ind,i1),fb_terms,d),1)' - D.*real(x(1:d,i1).^(gamma));
        end
    elseif isempty(fb_ind) % forward terms only
        for i1 = 1:size(x,2)
            f(1:d,i1) = S_sens*gf(:,i1) + S_fwd_*x(fwd_ind,i1) - D.*real(x(1:d,i1).^(gamma));
        end
    else % both feedback and forward terms exist
        for i1 = 1:size(x,2)
            f(1:d,i1) = S_sens*gf(:,i1) + sum(reshape(S_fb_.*x(fb_ind,i1),fb_terms,d),1)' + S_fwd_*x(fwd_ind,i1) - D.*real(x(1:d,i1).^(gamma));
        end
    end
    
    % GP dynamics
    f(gp_ind,:) = F*x(gp_ind,:);
    
end