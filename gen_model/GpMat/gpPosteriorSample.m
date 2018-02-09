function [fsamp, t] = gpPosteriorSample(kernType, numSamps, params, lims, observations ,seed, bw, sig_n, figNum)
  
% GPPOSTERIORSAMPLE Create a plot of samples from a posterior covariance.
% FORMAT
% DESC creates a plot of samples from a kernel with the given
% parameters and variance.
% ARG kernType : the type of kernel to sample from.
% ARG numSamps : the number of samples to take.
% ARG params : parameter vector for the kernel.
% ARG lims : limits of the x axis.
%
% COPYRIGHT : Neil D. Lawrence, 2008

% GP

global printDiagram

rand_obs = 0;
if nargin < 9
    figNum = [];
    if nargin < 8
      sig_n = 0;  
        if nargin < 7
          bw = false;
          if nargin < 6
            seed = [];
            if nargin < 5
                rand_obs = 1;
                if nargin < 4
                  lims = [-3 3];
                  if nargin < 3
                    params = [];
                    if nargin < 2
                      numSamps = 10;
                    end
                  end
                end
              end
            end
        end
    end
end
if ~isempty(seed)
  randn('seed', seed);
  rand('seed', seed);
end
t_star = linspace(lims(1), lims(2), lims(2))';

kern = kernCreate(t_star, kernType);
if isfield(kern, 'comp')
  for i=1:length(kern.comp)
    kern.comp{i}.transforms = [];
  end
end

if ~isempty(params)
  feval = str2func([kern.type 'KernExpandParam']);
  kern = feval(kern, params);
end
feval = str2func([kern.type 'KernExtractParam']);
[params, names] = feval(kern);
paramStr = [];
for i = 1:length(names)
  Name = names{i};
  Name(1) = upper(Name(1));
  ind = find(Name==' ');
  Name(ind+1) = upper(Name(ind+1));
  Name(ind) = '';
  paramStr = [paramStr Name num2str(params(i))];
  
end
paramStr(find(paramStr==46)) = 'p';
infoStr = ['Samples' num2str(numSamps) 'Seed' num2str(randn('seed'))];

% Covariance of the prior.
K_starStar = kernCompute(kern, t_star, t_star);

% Generate "training data" from a sine wave.
if rand_obs == 1
    t = rand(5, 1)*(lims(2)-lims(1))*0.6+lims(1)*0.6;
    f = sin(t);
else
    t = observations(:,1);
    f = observations(:,2);
end

% Compute kernel for training data.
K_starf = kernCompute(kern, t_star, t);
K_ff = kernCompute(kern, t);

%%%%% Added by WW %%%%%%
sig_n2I = diag(ones(size(K_ff,1),1)*(sig_n^2),0); % diagonal noise matrix sig_n^2I
K_ff = K_ff + sig_n2I; % fold in noise
%%%%%%%%%%%%%%%%%%%%%%%%

% Mean and covariance of posterior.
fbar = K_starf*pdinv(K_ff)*f;
Sigma = K_starStar - K_starf*pdinv(K_ff)*K_starf';

% Sample from the posterior.
fsamp = real(gsamp(fbar, Sigma, numSamps));
%{
% Plot and save
if isempty(figNum)
    figure
else
    figure(figNum);clf
end
linHand = plot(t_star, fsamp);
hold on
linHandPlot = plot(t, f, 'r.');
set(linHandPlot, 'markersize', 30)

zeroAxes(gca, 0.025, 18, 'times');
set(linHand, 'linewidth', 1)
app = '';
if bw
  set(linHand, 'color', [0 0 0])
  set(linHandPlot, 'color', [0 0 0])
  app = 'bw';
end

if iscell(kernType)
  KernType = [];
  for i = length(kernType):-1:1
    KernType = [kernType{i} KernType];
    KernType(1) = upper(KernType(1));
  end
else
  KernType(1) = upper(kernType(1));
end

if exist('printDiagram', 'var') & printDiagram
  printPlot(['gpPosteriorSample' KernType infoStr paramStr app], ...
            '../tex/diagrams', '../html')
end
%}

