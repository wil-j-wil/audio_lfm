% DEMINTERPOLATIONGP Demonstrate Gaussian processes for interpolation.
% FORMAT
% DESC runs a simple one-D Gaussian process displaying errorbars.
%
% SEEALSO : gpCreate, demRegressionGp
% 
% COPYRIGHT : Neil D. Lawrence, 2006, 2008

% GP

% Create data set
x = linspace(-1, 1, 9)';
trueKern = kernCreate(x, 'matern52');
K = kernCompute(trueKern, x);
% Sample some true function values.
yTrue = gsamp(zeros(size(x))', K, 1)';
%%
markerSize = 20;
markerWidth = 6;
markerType = 'b.';
lineWidth = 2;
% Create a test set
indTrain = [1 3 5 7 9]';
figNo = 1;
fillColor = [0.7 0.7 1];
for i = 0:1
  if i > 0
    yTrain = yTrue(indTrain);
    xTrain = x(indTrain);
    kern = kernCreate(x, 'matern52');
    % Change inverse variance (1/(lengthScale^2)))
    kern.inverseWidth = 5;
    
    xTest = linspace(-2, 2, 200)';
    xTestNew = 3000*(xTest+2);
    
    Kx = kernCompute(kern, xTest, xTrain);
    Ktrain = kernCompute(kern, xTrain, xTrain);
    
    yPred = Kx*pdinv(Ktrain)*yTrain;
    yVar = kernDiagCompute(kern, xTest) - sum(Kx*pdinv(Ktrain).*Kx, 2);
    ySd = sqrt(yVar);
  else
    p = [];
  end
  if i < length(indTrain)
    figure(16);clf
    if i>0
      fill([xTestNew; xTestNew(end:-1:1)], ...
           [yPred; yPred(end:-1:1)] ...
           + 2*[ySd; -ySd], ...
           fillColor,'EdgeColor',fillColor)
      hold on
      h = plot(linspace(1,12000,200), yPred, 'b-');
      %/~
      %      h = [h plot(xTest, yPred + 2*ySd, 'b--')];
      %      h = [h plot(xTest, yPred - 2*ySd, 'b--')];
      %~/
      set(h, 'linewidth', lineWidth)
    end
    p = plot(3000*(x(indTrain)+2), yTrue(indTrain), markerType);
    set(p, 'markersize', 20, 'linewidth', 1);
    ylim([-5 5])
    %{
    set(gca, 'xtick', [-2 -1 0 1 2]);
    set(gca, 'ytick', [-3 -2 -1 0 1 2 3]);
    set(gca, 'fontname', 'times', 'fontsize', 18, 'xlim', [-2 2], 'ylim', [-3 3])
    zeroAxes(gca);
    figNo = figNo + 1;
    %}
  end
end
