% MATERN_MODEL State space representation of GP with Matern kernel
%
% Syntax:
%   model = matern_model(theta,p,calc_grad)
% 
% In:
%        theta - GP parameters (length scale and magnitude) 
%        param - Static parameter structure with fields:
%        param.p - Integer part of the smoothness parameter nu = p + 1/2
%    calc_grad - Should the gradients be calculated
%     
% Out:
%        model - Structure containing the LTI model parameters and their
%                derivatives wrt given parameter vector theta               
% 
% Description:
%
%   Form a state space representation to Gaussian process with
%   Matern covariance function
%
%     C(t,t') = sigma^2 exp(-sqrt(2*(p+1)/2) tau/l) Gamma(p+1)/Gamma(2*p+1)
%             * sum_{i=0}^p [(p+i)!/(i!(p-i)!)(sqrt(8(p+1/2))tau/l)^(p-i)]
% 
%   where tau = |t-t'| and l is the length-scale (phi = sqrt(2)/2*l).
%

% Copyright (C) 2010-2012 Jouni Hartikainen, Simo Särkkä
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function model = matern_model(theta,param,calc_grad)
    
    if isnumeric(param)
        p = param;
    else
        p = param.p;
    end
    % Indexes to fixed parameters
    if isfield(param,'fixp_ind')
        fixp_ind = param.fixp_ind;            
    else
        fixp_ind = [];        
    end
    
    %
    if isfield(param,'fixp_val')
        fixp_val = param.fixp_val;
    end
    
    ntheta = 2;
    % Add the fixed parameter to the theta vector in right place
    k1 = 0; k2 = 0;
    thetaind = zeros(1,ntheta);
    for i = 1:ntheta
        if ismember(i,fixp_ind)
            k1 = k1 + 1;            
            theta = [theta(1:i-1) fixp_val(k1) theta(i:length(theta))];     
        else
            k2 = k2 + 1;
            thetaind(i) = k2; 
        end        
    end
    nthetao = ntheta-length(fixp_ind);
    
    if nargin < 3
       calc_grad = 0; 
    end

    model = struct;
    phi = sqrt(2)./2./theta(1);
    sigma = theta(2);
    la = 2*sqrt((p+0.5))*phi;
    c = sigma^2.*2*pi^(0.5)*exp(gammaln(p+1)-gammaln(p+0.5))*la.^(2*p+1);
   
    ppoly = pascal(p+2);
    ppoly = diag(flipud(ppoly))';
    ppoly = ppoly(2:end);
    
    lav = zeros(1,p+1);
    for i = 1:length(lav)
        lav(i) = la^i;
    end

    F = diag(ones(p,1),1);
    F(end,:) = fliplr(-lav.*ppoly); 
    
    L = zeros(p+1,1);
    L(end) = 1;
    H = zeros(1,p+1);
    H(1) = 1;

    Qc = L*c*L';
    
    % Prior mean and covariance
    M0 = zeros(size(F,1),1);
    %M0 = -5*ones(size(F,1),1); % WJW - Feb 2017
    try
        P0 = are(F',zeros(size(F)),Qc);
    catch
        theta
        error('No solution found for P0')
    end

    % Store to struct
    model.F  = F;
    model.G  = [];
    model.H  = H;
    model.Qc = Qc;
    model.M0 = M0;
    model.P0 = P0;
    
    % Gradients of model parameters 
    if calc_grad == 1
        % Gradient of F wrt theta
        lad = -sqrt(2*p+1)/(theta(1)^2);
        
        lavd = zeros(1,p+1);
        for i = 1:length(lavd)
            lavd(i) = i*la^(i-1)*lad;
        end
        DF = zeros(p+1,p+1,nthetao);

        if ~ismember(1,fixp_ind)
            DF(end,:,1) = fliplr(-lavd.*ppoly);
        end
        % Gradient of Qc wrt theta
        DQc = zeros(size(L,1),size(L,1),nthetao);

        if ~ismember(1,fixp_ind)
            DQc(:,:,thetaind(1)) = L*(2*p+1)*sigma^2.*2*pi^(0.5)*exp(gammaln(p+1)-gammaln(p+0.5))*la.^(2*p)*lad*L';
        end
        if ~ismember(2,fixp_ind)
            DQc(:,:,thetaind(2)) = L*2*sigma.*2*pi^(0.5)*exp(gammaln(p+1)-gammaln(p+0.5))*la.^(2*p+1)*L';
        end
        % Gradient of M0 wrt theta (is zero)
        DM0 = zeros(size(F,1),nthetao);

        % Gradient of P0 wrt theta
        DP0 = zeros(size(P0,1),size(P0,1),nthetao);

        for k = 1:2
            if ~ismember(k,fixp_ind)
                A = DF(:,:,thetaind(k))*P0 + P0*DF(:,:,thetaind(k))';
                DP0(:,:,thetaind(k)) = are(F',zeros(size(DF(:,:,thetaind(k)))),DQc(:,:,thetaind(k))+A);
            end
        end
        % Store to struct
        model.DF  = DF;
        model.DG  = [];
        model.DQc = DQc;
        model.DM0 = DM0;
        model.DP0 = DP0;
    end
    
end

