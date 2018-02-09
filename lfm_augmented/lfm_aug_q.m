% Edited version of TF_Q Dispersion matrix of transcription factor model
% from LFM toolbox
% 
% Syntax:
%   f = tf_q(x,t,param)
%
% In:
%   x     - State vector   (not used)
%   t     - Time instance  (not used)
%   param - Model parameters (see below for details)    
% Out:
%   Qc    - Dispersion matrix
%
% Description:
%   Dispersion matrix of transcriptional factor, described by the equations
%   
%   dx_i(t)/dt = B_i + \sum_{r=1}^R S_{i,r} g_i(u_r(t)) - D_i x(t),
% 
%   for outputs x_i(t), i=1,...,N and latent forces u_r(t), r=1,...,R,
%   which modelled with LTI SDE models specified by the user in param
%   structure. The non-linear function g_i(.) is also specified in param.
%
%   Currently the output model has no stochasticity, so the only non-zero
%   components of the dispersion matrix belong to force components.
% 
% Initially based on code by Jouni Hartikainen, Simo Särkkä
%
% Modified for audio LFM work by Will Wilkinson - Jan 2018

function [Qc] = lfm_aug_q(x,t,param)
    
    % Output q
    d = param.d;
    Qc = zeros(d);
    
    % LFM q
    Qgp = sqrt(param.Qgp);
    Qc = blkdiag(Qc,Qgp);
    
    % Feeback & Forward q
    %fb_terms = param.fb_terms;
    %fwd_terms = param.fwd_terms;
    %Qf = zeros(d*fb_terms+fwd_terms);
    %Qc = blkdiag(Qc,Qf);

end