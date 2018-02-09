function GPsamp = stateSpaceGpSample(lengthscale,sigma,time_steps)

%% Sample from state-space GP given hyperparameters

%% Parameters
%lengthscale = 5; % resultsMag.theta_opt(end-1) = 2.7803
%sigma = 0.1; % resultsMag.theta_opt(end) = 2.1747 (0.05 seems more suitable (for dt=0.1). Why?)
blockSize = 512;
numBlocks = ceil(time_steps/blockSize);

dt = 0.1;% * (342 / 175906);

%% Matern GP type
p = 2;

%% Model calcs.
%tic
phi = sqrt(2)./2./lengthscale;
la = 2*sqrt((p+0.5))*phi; % nu = p + 0.5
% see eqn (3.12) on p51 of Hartakeinen's thesis (2013)
% "spectral density of white noise process"
% gammaln is the log of the gamma function, hence the exponent
c = sigma^2.*2*pi^(0.5)*exp(gammaln(p+1)-gammaln(p+0.5))*la.^(2*p+1); % q

ppoly = pascal(p+2);
ppoly = diag(flipud(ppoly))';
ppoly = ppoly(2:end);

lav = zeros(1,p+1);
for i = 1:length(lav)
    lav(i) = la^i;
end

% F and L are the matrices used for companion form of the LTI SDE
F_lf = diag(ones(p,1),1);
F_lf(end,:) = fliplr(-lav.*ppoly);

L = zeros(p+1,1);
L(end) = 1;
H_lf = zeros(1,p+1);
H_lf(1) = 1;

Qc_lf = L*c*L'; % (3.7) on p50 of Hartikeinen Thesis (2013)

%% Discretise ([Alf,Qlf] = lti_disc(F_lf,[],Qc_lf,dt);)
% Closed form integration of transition matrix
Alf = expm(F_lf*dt);

% Closed form integration of covariance by matrix fraction decomposition
n_lf = size(F_lf,1);
L_lf = eye(size(F_lf,1));
Phi = [F_lf L_lf*Qc_lf*L_lf'; zeros(n_lf,n_lf) -F_lf'];
AB = expm(Phi*dt)*[zeros(n_lf,n_lf);eye(n_lf)];
Qlf = AB(1:n_lf,:)/AB((n_lf+1):(2*n_lf),:);

%% Draw sample from the GP

%% Setup
%figure
x0 = zeros(size(Alf,1),1); % Inital conditions
x = x0;
t = 1:blockSize;
cholQlf = chol(Qlf)';
cholQlfvec = cholQlf*ones(size(cholQlf,1),1);
X = zeros(3,blockSize);
Xlf = zeros(1,blockSize);
%toc

%% Block by block sample generation
GPsamp = zeros(blockSize*numBlocks);
gind = 1:blockSize;
for n = 1:numBlocks
    %tic
    for i = 1:blockSize;
        x = Alf*x + cholQlfvec*randn(1); % Share the same noise sample across dimensions
        X(:,i) = x;
    end
    Xlf = H_lf*X;
    GPsamp(gind) = Xlf';
    gind = gind + blockSize;
    %toc
    %plot(t + (n-1)*blockSize, Xlf)
    %hold on
end
%title('State-Space (Kalman Filtering) Sample of Latent Force')

GPsamp = GPsamp(1:time_steps);

end