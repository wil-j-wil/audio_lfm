function [tNMF_w,tNMF_h] = run_tNMF(ATrain,K)

% A function to run all the stages of Richard Turner's temporal NMF (tNMF)
% code

% K = # latent activations / features / functions / whatever you want to call them

[T,D] = size(ATrain);

%% Train the model
HInit = exp(randn(T,K));
ks = ceil(T*rand(K,1));
WInit = ATrain(ks,:);
vary = zeros(T,D);

% initialise with regular NMF
Opts.restarts = 50;
Opts.numIts = 5000;
%[WEst1,HEst1,info1] = nmf(ATrain,WInit,HInit,[],[],[],vary,Opts);
[WEst1,HEst1,~] = nmf_fp(ATrain,WInit,HInit,vary,Opts);
%[HEst1,WEst1] = nnmf(ATrain,K);

% order according to slowness (mean square derivative)
fastness = mean(diff(HEst1).^2)./var(HEst1);
[~,ind] = sort(fastness,'descend');
HEst1 = HEst1(:,ind); 
WEst1 = WEst1(ind,:);

%% LOG-NORMAL SE GAUSSIAN PROCESS TEMPORAL NMF
% initialise

lenx = zeros(K,1);
mux = zeros(K,1);
varx = zeros(K,1);

for k=1:K
  % threshold
  logHthresh = log(HEst1(:,k)+1e-8);
  
  % filter
  filt = exp(-1/2*([-100:100].^2)/(1^2));
  filt = filt/sum(filt);
  logHsm = conv(logHthresh,filt,'same');
  
  % fit GP
  mux(k) = mean(logHsm);
  [lenx(k),varx(k),~] = trainSEGP_RS(logHsm-mux(k));
  
end

disp(['length scale parameters'])
lenx
Opts.numIts = 200;
[HEst2,~] = tnmf_inf(ATrain,WEst1,HEst1+1e-8,lenx,mux,varx,vary,Opts);
Opts.numIts = 2000;
[WEst2,HEst2,~] = tnmf(ATrain,WEst1,HEst2,lenx,mux,varx,vary,Opts);
% confusingly, Turner uses W and H the opposite way round to Matlab's nnmf
% function, so I do a switch here:
tNMF_w = HEst2;
tNMF_h = WEst2;
end