function [envSynth,latentForces,sigSynth,carriers,latentForcesDS] = reconstruct_sig(lfmData,lfSynth,RemoveNoiseFloor,plotResults,newLength,whichForce,residual,initCond)

% takes the optimised model, stored in 'lfmData', and uses the latent
% forces stored at 'lfmData.filterFinal.uAdjusted' or some user provided
% forces to generate amplitude envelopes.

% Will Wilkinson - Jan 2018

if nargin < 2
    lfSynth = [];
end
if nargin < 3 || isempty(RemoveNoiseFloor)
    RemoveNoiseFloor = 0; % set to 1 to remove noise floor
end
if nargin < 4 || isempty(plotResults)
    plotResults = 0; % set to N to plot results on figure N
end
if nargin < 5 || isempty(newLength)
    newLength = 0;
end
if nargin < 6
    whichForce = [];
end
if nargin < 7
    residual = [];
end
if nargin < 8 || isempty(initCond)
    initCond = lfmData.filterFinal.At;
end

if newLength == 0
    lf_targetSize = size(lfmData.filterFinal.Y,2);
    out_targetSize = length(lfmData.demod_struct.y);
else
    lf_targetSize = max(size(lfSynth));
    out_targetSize = max(size(lfSynth))*lfmData.lfmOpt.DS;
end

if isempty(lfSynth) == 0
    forces = log(exp(lfSynth)-1);
    if isempty(whichForce)
        forcesToExclude = [];
    else
        forcesToExclude = setdiff(1:lfmData.lfmOpt.nlf,whichForce);
    end
    for j=intersect(1:lfmData.lfmOpt.nlf,forcesToExclude)
        forces(j,:) = -100; % i.e. pretty much zero when transformed
    end
    lf_synth = makeThisThatLong(forces,lf_targetSize);
else
    forces = lfmData.filterFinal.uAdjusted;
    if isempty(whichForce)
        forcesToExclude = [];
    else
        forcesToExclude = setdiff(1:lfmData.lfmOpt.nlf,whichForce);
    end
    for j=intersect(1:lfmData.lfmOpt.nlf,forcesToExclude)
        forces(j,:) = -100; % i.e. pretty much zero when transformed
    end
    lf_synth = makeThisThatLong(forces,lf_targetSize);
end

mods = lfmData.demod_struct.A; % modulators

allU = [];
allUDS = [];

ssOut = runSS(lfmData.filterFinal,[],lf_synth,[],0,initCond);
U = log(1+exp(lf_synth));
ssOutUp = makeThisThatLong(ssOut,out_targetSize);
UUp = makeThisThatLong(U,out_targetSize);
UUp(find(UUp<0)) = 0; %#ok<FNDSB>

envSynth = ssOutUp';
allUDS = [allUDS; log(1+exp(lf_synth))];
allU = [allU; UUp];

if ~isempty(residual)
    ssmax = max(ssOut)';
    residual = bsxfun(@times,residual,ssmax);
    residual = makeThisThatLong(residual,out_targetSize);
    envSynth = envSynth + residual;
end

if RemoveNoiseFloor > 0
    perc = zeros(1,size(envSynth,2));
    percm = zeros(1,size(envSynth,2));
    for j=1:size(envSynth,2)
        perc(j) = prctile(envSynth(:,j),RemoveNoiseFloor); % compare ?th percentile
        percm(j) = prctile(mods(:,j),RemoveNoiseFloor);
    end
    percDiff = max(perc - percm,0);
    envSynth = bsxfun(@plus,envSynth,-percDiff);
end
envSynth(find(envSynth<0))=0; %#ok<FNDSB>
latentForces = allU;
latentForcesDS = allUDS;

if newLength > 0
    envSynth = makeThisThatLong(envSynth,newLength);
    latentForces = makeThisThatLong(latentForces,newLength);
end

envSynth = envSynth / lfmData.env_scale;
envSynth(:,lfmData.I) = envSynth;

if plotResults > 1
    figure(plotResults);clf
    subplot(311)
    plot(lfmData.demod_struct.A)
    title('Actual Amplitude Envelopes')
    subplot(312)
    plot(allU','LineWidth',2)
    title('Latent Forces')
    subplot(313)
    plot(envSynth)
    title('Reconstructed Amplitude Envelopes (Upsampled)')
    xlabel('time (samples)')
end

carrier_synth = lfmData.demod_struct.carrier_synth;
if newLength == 0
    carriers = carrier_synth(1:size(envSynth,1),:);
else
    carriers = zeros(size(envSynth));
    for i=1:size(envSynth,1)
        carriers(i,:) = carrier_synth(mod(i-1,size(carrier_synth,1))+1,:);
    end
end

sigSynth = sum(carriers'.*envSynth');

end