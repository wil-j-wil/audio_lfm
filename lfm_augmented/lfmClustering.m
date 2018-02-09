function [cidx,clustTreeCos,distances,kmeansDist,numClusters] = lfmClustering(envelopes,numClusters,shiftMax,clustMinSize,clustMaxSize,numIters)

% this function is used to cluster envelopes based on cosine similarity -
% used to set initial parameters prior to LFM optimisation.

% Will Wilkinson - Jan 2018

numOut = size(envelopes,1);

if nargin < 4
    clustMinSize = 1;
end
if nargin < 5
    clustMaxSize = numOut;
end
if nargin < 6
    numIters = 100;
end

cosD = ones(numOut);
Dind = zeros(numOut);
for i=1:numOut
    for j=i+1:numOut
        for k=1:2*shiftMax
            envelope1 = envelopes(i,shiftMax+1:end-shiftMax);
            envelope2 = envelopes(j,k:end-(2*shiftMax)-1+k);
            cosDist = pdist([envelope1;envelope2],'cosine');
            if cosDist < cosD(i,j)
                cosD(i,j) = cosDist;
                Dind(i,j) = k-shiftMax;
            end
        end
    end
end

envelopes_ = envelopes;
D = size(envelopes,1);
for i=1:D
    if Dind(1,i) > 0
        envelopes_(i,:) = [envelopes(i,Dind(1,i)+1:end) zeros(1,Dind(1,i))];
    elseif Dind(1,i) < 0
        envelopes_(i,:) = [zeros(1,-Dind(1,i)) envelopes(i,1:end+Dind(1,i))];
    end
end

% scale the magnitudes
%adjustFactor = 0.1 ./ max(envelopes_,[],2);
%envelopes_ = bsxfun(@times,envelopes_,adjustFactor);
[cidx,clustTreeCos,distances,kmeansDist] = partialCluster(envelopes_,numClusters,numIters,'cosine');

numels = zeros(numClusters,1);
for i=1:numClusters
    numels(i) = length(find(cidx==i));
end

if clustMaxSize*numClusters >= numOut;
    while max(numels) > clustMaxSize
        [~,clusterRank] = sort(numels,'descend');
        els = find(cidx==clusterRank(1));
        compare = 1:numClusters;
        compare(find(numels>=clustMaxSize)) = [];
        [D,Ir] = min(kmeansDist(els,compare),[],2);
        [~,Ic] = min(D);
        toMove = els(Ic);
        moveTo = compare(Ir(Ic));
        cidx(toMove) = moveTo;
        for i=1:numClusters
            numels(i) = length(find(cidx==i));
        end
    end

    stop = 0;
    while min(numels) < clustMinSize && stop < 1
        [~,clusterRank] = sort(numels,'ascend');
        els = find(cidx==clusterRank(1));
        compare = 1:numClusters;
        compare([clusterRank(1); find(numels>=clustMaxSize)]) = [];
        if isempty(compare)
            stop = 1;
        else
            [D,Ir] = min(kmeansDist(els,compare),[],2);
            [~,Ic] = min(D);
            toMove = els(Ic);
            moveTo = compare(Ir(Ic));
            oldClust = cidx(toMove);
            cidx(toMove) = moveTo;
            find(cidx==oldClust);
            if isempty(find(cidx==oldClust))
                cidx(find(cidx>oldClust)) = cidx(find(cidx>oldClust)) - 1;
                numClusters = numClusters - 1;
                fprintf('number of clusters was reduced to %d \n',numClusters)
            end
            numels = zeros(numClusters,1);
            for i=1:numClusters
                numels(i) = length(find(cidx==i));
            end
        end
    end
else
    disp('clustering configuration not possible')
end

end