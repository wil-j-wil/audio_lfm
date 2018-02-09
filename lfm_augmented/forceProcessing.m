function lfmDataOut = forceProcessing(lfmDataIn,plotNum)

% passing our Gaussian predictive dist. through the nonlinearity gives our
% predictions an upwards bias - so here we adjust down based on which
% percentile gives us the smallest reconstruction error

% Will Wilkinson - Jan 2018

if nargin < 2
    plotNum = 0;
end

newU = lfmDataIn.filterFinal.ForceSMAP;
G = log(1+exp(newU));

ssOut_ = runSS(lfmDataIn.filterFinal,0,[],[],0);

inSig = lfmDataIn.filterFinal.inSig;
sigMax = max(inSig)/max(max(inSig));
for l=1:lfmDataIn.lfmOpt.nlf-1
    sigMax = [sigMax; max(inSig)/max(max(inSig))]; %#ok<AGROW>
end

actEnv = lfmDataIn.filterFinal.inSig;

uVar = zeros([size(lfmDataIn.filterFinal.uVarF) 2]);
uVar(:,:,1) = lfmDataIn.filterFinal.uVarS;
uVar(:,:,2) = lfmDataIn.filterFinal.uVarF;
strength_rng = 0:0.5:10;
skew_rng = [0.5 0.75 1:0.5:20];
error_store = 999*ones(length(strength_rng),length(skew_rng),2);
for k=1:2
    i=1;
    for adjustStrength=strength_rng % strength
        j=1;
        for skew=skew_rng
            adjustedU = newU - adjustStrength*sqrt(abs(uVar(:,:,k)));
            lowLim = adjustedU;
            Glow = log(1+exp(lowLim));

            G_ = G/max(max(G));
            Glow_ = Glow/max(max(G));

            %adjustedG_ = G_ - (1-G_).^skew.*(G_-Glow_);
            adjustedG_ = G_ - ((1-G_).^skew + (1-sigMax).^skew)./2 .*(G_-Glow_);
            adjustedG = adjustedG_ * max(max(G));
            adjustedU_ = log(exp(adjustedG)-1);

            ssOut = runSS(lfmDataIn.filterFinal,0,adjustedU_,[],0);
            errLFM_adjusted = sum(sum((ssOut-actEnv).^2));
            error_store(i,j,k) = errLFM_adjusted;
            j=j+1;
        end
        i=i+1;
    end
end

[~,min_error_k] = min(min(min(error_store)));
[~,min_error_j] = min(min(error_store(:,:,min_error_k)));
[~,min_error_i] = min(error_store(:,min_error_j,min_error_k));

%error_store(min_error_i,min_error_j,min_error_k)
ad_strg = strength_rng(min_error_i);
skw = skew_rng(min_error_j);
uVar = uVar(:,:,min_error_k);

adjustedU = newU - ad_strg*sqrt(abs(uVar));
lowLim = adjustedU;
Glow = log(1+exp(lowLim));

G_ = G/max(max(G));
Glow_ = Glow/max(max(G));

%adjustedG_ = G_ - (1-G_).^skw.*(G_-Glow_);
adjustedG_ = G_ - ((1-G_).^skw + (1-sigMax).^skw)./2 .*(G_-Glow_);
adjustedG = adjustedG_ * max(max(G));
adjustedU_ = log(exp(adjustedG)-1);

ssOut = runSS(lfmDataIn.filterFinal,0,adjustedU_,[],0);
errLFM = sum(sum((ssOut_-actEnv).^2));
errLFM_adjusted = sum(sum((ssOut-actEnv).^2));
fprintf('LFM error: %s\n',errLFM)
fprintf('Adjusted LFM error: %s\n',errLFM_adjusted)

t = 1:size(newU,2);

color1 = [0 0 1];
color1_ = [0.7 0.7 1];

uncertaintyUp = newU+1.96*sqrt(abs(uVar));
uncertaintyLow = newU-1.96*sqrt(abs(uVar));

if plotNum > 0
    figure(plotNum);clf
    subplot(211)
    h1=fill([t fliplr(t)], [uncertaintyUp fliplr(uncertaintyLow)], color1_, 'edgecolor',color1_,'DisplayName','Uncertainty');
    hold on
    h2=plot(newU','color',color1,'DisplayName','Mean');
    %plot(adjustedU','r')
    h3=plot(adjustedU_','r','DisplayName','Adjusted Mean');
    legend([h1(1) h2(1) h3(1)])
    title('The Predicted Latent Force, Mean Adjusted')
    subplot(212)
    h1=fill([t fliplr(t)], [log(1+exp(uncertaintyUp)) fliplr(log(1+exp(uncertaintyLow)))], color1_, 'edgecolor',color1_,'DisplayName','Uncertainty');
    hold on
    h2=plot(G_'*max(max(G)),'color',color1,'DisplayName','Mean');
    %plot(Glow_'*max(max(G)),'r')
    h3=plot(adjustedG_'*max(max(G)),'r','DisplayName','Adjusted Mean');
    legend([h1(1) h2(1) h3(1)])
    title('The Postitive Predicted Latent Force, Mean Adjusted')
end

lfmDataOut = lfmDataIn;
lfmDataOut.filterFinal.uVar = uVar;
lfmDataOut.filterFinal.uAdjusted = adjustedU_;

end