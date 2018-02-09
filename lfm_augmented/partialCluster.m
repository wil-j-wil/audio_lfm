function [cidx,clustTreeCos,distances,kmeansDist] = partialCluster(partials,k,reps,distType)

% kmeans clustering function

% Will Wilkinson - Jan 2018

if nargin < 3
    reps = 10;
end
if nargin < 4
    distType = 'cosine';
end

%% hierarchical
distances = pdist(partials,distType);
clustTreeCos = linkage(distances,'average');

% ^ observe the dendrogram and pick a value for k

%% k-means

[cidx,~,~,kmeansDist] = kmeans(partials,k,'dist',distType,'replicates',reps,'display','off'); % cidx is our clustering index


