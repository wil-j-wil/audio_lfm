% Dynamic model for moments with Gaussian filters (square root version)
% 
% Syntax:
%  dmAC = dmAC_sqrt(vmAC,param)
%
% In:
%   vmAC  - m, A (sqrt of P, P = AA') and cross covariance C in a single vector
%   param - Parameter vector (see below for details)
%      
% Out:
%   dmAC - Time derivative of M, A and C
%
% Initially based on code by Jouni Hartikainen, Simo Särkkä
%
% Modified for audio LFM work by Will Wilkinson - Jan 2018


function dmAC = lfm_aug_dmAC_sqrt(vmAC,param)    
    XI        = param{1};     % Unscaled sigma points
    WM        = param{2};     % Mean weights of sigma points
    WC        = param{3};      % Covariance weights of sigma points
    f_func    = param{4};     % Dynamic model function
    f_param   = param{5};     % Static parameters of f
    Q         = param{6};     % Diffusion matrix 
    t         = param{end};   % Time instance

    n = size(XI,1);
    n_l = f_param.d*f_param.fb_terms + f_param.nlf*f_param.fwd_terms;

    [mAC,m,A,C] = vec2mAC(vmAC,n+n_l);
     
    XI = [XI; zeros(n_l,size(XI,2))];
    % Multiply with Cholesky and add mean
    %X = A*XI + repmat(m,1,size(X,2));
    % Faster way:
    X = bsxfun(@plus,A*XI,m);
    
    % Evaluate the function at the sigma points
    Z = feval(f_func,X,t,f_param);
    Z_l = Z(n+1:n+n_l,:);
    Z = Z(1:n,:);
    
    % Evaluate Q if it is not fixed
    if length(param) > 7
        q_func  = param{7};
        q_param = param{8};

        Qc = feval(q_func,X,t,q_param);
        Q = zeros(n,n);
        
        % Approximate the expectation E[L(x,t) Qc L(x,t)']
        for i = 1:size(Z,2)
            Q = Q + WC(i).*Qc(:,:,i)*Qc(:,:,i)';
        end
    end
    
    % Approximate E[f(x,t)]
    Ef = zeros(n,1);
    for i=1:size(Z,2)
       Ef = Ef + WM(i).*Z(:,i); 
    end
    
    % ... and E[(x-m)f(x,t)^T]
    D  = zeros(n,size(Ef,1));
    for i=1:size(Z,2)
        D  = D + WC(i) * (X(1:n,i) - m(1:n)) * (Z(:,i) - Ef)';
    end
   
    dm = Ef;
    dm_l = zeros(n_l,1); % zeros for - change to bucket brigade
    dm = [dm; dm_l];

    B = D+D'+Q;
    A = A(1:n,1:n);
    B = (A\B)/A';
    
    B = tril(B) - 1/2*diag(diag(B));
    dA = A*B;

    % dC/dt = C P^-1 E[(x-m)f(x,t)']
    dC = C(1:n,1:n)*(A'\(A\D));
    dC_l = zeros(n_l); % zeros for - change to bucket brigade
    dC = blkdiag(dC,dC_l);
    
    dA_l = zeros(n_l); % zeros for - change to bucket brigade
    dA = blkdiag(dA,dA_l);
    dA = dA(logical(tril(ones(n+n_l))));
    dmAC = [dm(:); dA(:); dC(:)];
    
