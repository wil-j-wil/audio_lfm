function ssOut = runSS(results,displayParam,U,T,plotFig,initCond)

% Run the state space model

% Will Wilkinson - Jan 2018

N = size(results.inSig,1);
nlf = results.nlf;
sparsity_fb = results.sparsity_fb;
sparsity_fwd = results.sparsity_fwd;
fb_terms = length(sparsity_fb);
fwd_terms = length(sparsity_fwd);
%gamma = results.gamma;
dt = results.dt;

if nargin < 2 || isempty(displayParam)
    displayParam = 0;
end
if nargin < 3 || isempty(U)
    U = results.ForceSMAP;
end
if nargin < 4 || isempty(T)
    T = size(U,2);
end
if nargin < 5 || isempty(plotFig)
    plotFig = 99;
end
if nargin < 6 || isempty(initCond)
    initCond = results.At;
end

if plotFig > 0
    figure(plotFig);clf
    subplot(311)
    plot(results.inSig')
    %actualAxis=[xlim ylim];
    title('Actual Data')
    subplot(312)
    plot(log(1+exp(U(1,:)))','k--','LineWidth',2)
    hold on
    plot(log(1+exp(U(2:end,:)))','--','LineWidth',2)
    title('Latent Function(s)')
end

%A = results.At; % Initial state
A = initCond; % Initial state
B = results.Bt; % Offset

[~, gamma, D_, ~, ~, S_sens, S_fb, S_fwd] = unpackTheta(results);

D = eye(N,N);
D(eye(N)==1) = -D_; % Damping

if displayParam == 1 %#ok<*NOPRT>
    % display parameters:
    D
    S_sens
    S_fb
    S_fwd
end

G = log(1+exp(U));

%G(84:end) = 0;

% run the state space model
y = zeros(N,T);    % Preallocate output signal for t=1:T

% Perform the system simulation:
x = A;                                  % Set initial state
for t=1:T                               % Iterate through time
    y(:,t) = x;% .* adjustFactor';      % Output for time t
    SG = S_sens*G(:,t);
    for k=1:min(fb_terms,t-1)
        SG = SG + S_fb(:,k).*y(:,t-k); % feedback
    end
    for k=1:min(fwd_terms,t-1)
        SG = SG + S_fwd(:,k:fwd_terms:fwd_terms*nlf)*G(:,t-k); % forward
    end
    %t
    %D*(x.^gamma)
    %SG
    %B
    xdot = D*(x.^gamma) + SG + B;       % calculate x'(t)
    x = x + (xdot * dt);                % calculate x(t+1) = forward one time step
    x(find(x<0)) = 0; %#ok<*FNDSB>
end

% adjust the magnitudes (this is allowed since it's just an additional scalar operation)
%maxy = max(y,[],2);
%maxY = max(results.inSig,[],2);
%adjustFactor = maxY ./ maxy;
%y = bsxfun(@times,y,adjustFactor);

if exist('results.scaleFactor','var') == 1
    ssOut=y/results.scaleFactor;
else
    ssOut=y;
end

if plotFig > 0
    figure(plotFig);
    subplot(313)
    plot(ssOut')
    title('Synthesised Output')
    %axis(actualAxis)
    drawnow()
end

end