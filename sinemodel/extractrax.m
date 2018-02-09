function [trax,mags,dbgs] = extractrax(s,h,verbose)
% [T,M] = extractrax(S,H,THR,VERB)  Extract tracks from a 2-d t-f style array
%       S is a matrix of values, which is searched from left to right 
%       looking for local maxima in each column, then tracking their 
%       evolution.  T is returned with one row for each track that 
%       is found, with nonzero values indicating the continuously-
%       valued  row (i.e. frequency, but starting from 1.0) of the 
%       track in the original matrix.  
%       H is a threshold; only peaks > H*maxx(S) at startup 
%       are tracked. H defaults to 0.01.
%       Rows of T are sorted by column (i.e. start time), then by freq, 
%       of their first value. VERB=1 for progress printouts.
%	M returns the interpolated peak magnitude for each track point.
% 1998may02 dpwe@icsi.berkeley.edu $Header: $ how many times?

% Performance constants:

if nargin < 2
  % Default proportion-of-max startup threshold cutoff
  h = 0.005;
end
if nargin < 3
  verbose = 0;
end

  
% Ratio of peak-height-to-start vs peak-height-to-continue (hysteresis)
startupmargin = 5.0;

% Maximum col-to-col tolerable peak movement (in bins)
maxdfdt = 3.0;   % used to be 1.0

% How many steps before say a peak is finished
deadsteps = 3;

% However, energy in the channel below this factor times the last 
% seen peak means it really has gone
lowEthresh = 0.2;

% Prune away tracks that end up with fewer points than this
minNpts = 10;

% Maximum mag-increase before new track (dB/step)
maxdBincr = 20.0;

% Asymptote of adaptive magnitude-increase threshold level (dB/step)
finaldBincr = 1.5; % was 6.0

% Adptv mag thresh stays this far above allowed positive slopes (must be > 1.0)
dbIncrSafeRatio = 1.5;

% Decay time for adptv mag incr threshold (steps)
dbIncrTsteps = 20;

%%%

% Figure the absolute threshold
startupThresh = h*max(max(s));

[nr,nc] = size(s);

%% % Track all the local maxes we looked at
%% d = zeros(size(s));

% 'frequencies' (i.e. bin indexes) of current track's most recent peaks
Fs = [];
% 'energies' (i.e. values) of current track's most recent peaks
Es = [];
% Start time of each track (to calc decaying mag thresh)
Ss = [];
% 'times' (column indices) of current track's most recent peaks 
% this can be less than the previous column because of permitted bridging
Ts = [];
% row numbers within the output "trax" array for the data records 
% of each current track
Rs = [];
% Current magnitude threshold limit.  This adapts
Ms = [];

blockrows = 50;
tblock = NaN*ones(blockrows, nc);
ntrax = 0;
% initialize the data array
trax = tblock;
% Also follow our extracted magnitudes, for yuks
mags = tblock;
% Debug parameters the same size
dbgs = tblock;

txsize = size(trax, 1);

% Do the columnwise search
for c = 1:nc
  col = s(:,c)';

  % find locations and heights of best local maxes
  % (return peaks even startupmargin smaller that startupThresh, 
  %  since that's the acceptable continuation height; filter later)
% original: fit peaks in linear domain
%  [maxFs, maxEs] = quadmaxloc(col, startupThresh/startupmargin);
% average |sgram| dB difference, weighted by resynth |sgram| = 0.3585 dB
% new: fit peaks in dB domain:
  [maxFs, maxEs] = quadmaxloc(log(col), log(startupThresh/startupmargin));
  maxEs = exp(maxEs);
% average |sgram| dB difference, weighted by resynth |sgram| = 0.0868 dB
% - big improvement!

  % Match up to existing tracks
  numcurrent = size(Fs, 2);
  % Setup mask to indicate which tracks should be terminated this time
  allowCont = ones(1,numcurrent);
  for i = 1:numcurrent
   % fprintf(1, 'col %d: extending current#%d (r=%d) f=%f\n', c, i, Rs(i), Fs(i));
    f = Fs(i);
    foundcont = 0;
      % Unfortunately, when maxFs is empty, this test crashes, so have to 
      % do it as 2 nested tests.
      %if size(maxFs,2) > 0 & abs(f - maxFs(bestix)) < maxdfdt
    if size(maxFs,2) > 0
      dfs = abs(maxFs - f);
      bestix = find(dfs == min(dfs));
      if size(maxFs,2) > 0 & abs(f - maxFs(bestix)) < maxdfdt
        % Found a continuation for this current peak
        % Check that it doesn't fail the magnitude increase check
%        magthr = finaldBincr + (maxdBincr-finaldBincr)*exp(-(c - Ss(i))/dbIncrTsteps);
	magthr = Ms(i);
	magdBStep = 20*log10(maxEs(bestix)/Es(i));
	if magdBStep > Ms(i)
	  % Step was too big: Kill this one & leave the continuation unused
	  % so that a new track is created with it
          if verbose ~= 0
            fprintf(1, ['col %d (l=%d) tk %d: large step of %f at ' ...
            '%f, breaking track\n'], c, c-Ss(i), i, magdBStep, f)
          end
	  allowCont(i) = 0;
	else
          Fs(i) = maxFs(bestix);
          Es(i) = maxEs(bestix);
          Ts(i) = c;
          % delete the taken peak fro maxFs & maxEs
          mask = [1:size(maxFs,2)] ~= bestix;
          maxFs = maxFs(find(mask));
          maxEs = maxEs(find(mask));
          % Store 'freq' in output array
          trax(Rs(i), c) = Fs(i);
          mags(Rs(i), c) = Es(i);
	  dbgs(Rs(i), c) = Ms(i);
          foundcont = 1;
	  % Update the mag threshold
	  Ms(i) = max(dbIncrSafeRatio*magdBStep, ...
		      Ms(i) - (Ms(i)-finaldBincr)*(1-exp(-1/dbIncrTsteps)));
        end
      end
    end
    if foundcont == 0
      % No acceptable continuation found - maybe kill this track
      if (((c - Ts(i)) > deadsteps) | (col(round(Fs(i))) < lowEthresh*Es(i)))
        % Mark this track to be removed from current structures
        allowCont(i) = 0;
      end
      % otherwise, let it hang on for a bit
      % Set a value for the frequency anyway?  Else could leave gaps
      % in frq to indicate untracked parts...
      %trax(Rs(i), c) = Fs(i);
      %mags(Rs(i), c) = Es(i);
    end
  end % of loop over currently active tracks

  % Remove the tracks that have ended from the current records
  if min(allowCont) == 0 
    Fs = Fs(find(allowCont));
    Es = Es(find(allowCont));
    Ts = Ts(find(allowCont));
    Rs = Rs(find(allowCont));
    Ss = Ss(find(allowCont));
    Ms = Ms(find(allowCont));
  end

  % Create new tracks for all remaining peaks
  for i = 1:(size(maxFs,2))
    if maxEs(i) > startupThresh
      if verbose ~= 0
        fprintf(1, 'Creating track at col %d bin %f\n', c, maxFs(i));
      end
      ntrax = ntrax + 1;
      tkrow = ntrax;
      if ntrax > txsize
        % Need to expand the pre-allocated trax buffer
	trax = [trax; tblock];
	mags = [mags; tblock];
	dbgs = [dbgs; tblock];
	txsize = size(trax,1);
        if verbose ~= 0
          fprintf(1, 'trax extended to %d rows\n', txsize);
        end
      end
      % trax = [trax; NaN*ones(1,nc)];
      % mags = [mags; NaN*ones(1,nc)];
      % dbgs = [dbgs; NaN*ones(1,nc)];
      % tkrow = size(trax,1);
      Rs = [Rs, tkrow];
      Fs = [Fs, maxFs(i)];
      Es = [Es, maxEs(i)];
      Ts = [Ts, c];
      Ss = [Ss, c];
      Ms = [Ms, maxdBincr];
      j = size(Rs,2);
      trax(Rs(j), c) = Fs(j);
      mags(Rs(j), c) = Es(j);
      dbgs(Rs(j), c) = Ms(j);
    end
  end

end

% Prune away tracks with too few pts (including pre-allocated empty ones!)
keeprows = sum(~isnan(trax)') >= minNpts;
trax = trax(find(keeprows),:);
mags = mags(find(keeprows),:);
dbgs = dbgs(find(keeprows),:);

% Offset freqs to start from zero, not 1 (2002-03-04)
trax = trax - 1;

% foreach col
%   find the local maxes
%   foreach local max
%     find the 'exact' location
%   foreach current track
%     if (nearest new peak)
%       extend track
%       remove that peak
%     elseif bridging less than max and actual val not too small
%       guess-extend-track 
%     else
%       remove track from active
%   foreach remaining new peak
%     create new active track
% done


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xmax,ymax] = quadmaxloc(v,th)
% [X,Y] = quadmaxloc(V,T)  Quadratic-interpolated index and values for loc maxs
%      V is a uniformly-sampled vector.  Find all local maxes above T 
%      (absolute), then do a quadratic fit to interpolate the location and 
%      height of the maxima.   Return these as correspoding elements of X
%      and Y.
% 1998may02 dpwe@icsi.berkeley.edu $Header: $ again?

if (nargin < 2) 
  th = 0;
end

if (size(v,1) > 1) 
  error('v must be a row vector');
end

nr = size(v,2);

% catch for error case
if max(v) == -Inf
  % whole vector is -Inf, i.e. frame was digital zero
  xmax = [];
  ymax = [];
else
  % filter for local maxima; ensure edges don't win
  gtl = (v > [v(1), v(1:(nr-1))]);
  % allow greater-than-or-equal to catch plateaux
  gtu = (v >= [v(2:nr), 1+v(nr)]);
  vmax = v .* (v > th) .* gtl .* gtu;
  maxixs = find(vmax);
  % Interpolate the max pos's
  xmax = zeros(size(maxixs));
  ymax = zeros(size(maxixs));
  for i = 1:size(maxixs,2)
    % Solve quadratic fit to 3 pts (as y = ax(x-b) with 0,0 as col(rmax-1))
    rmax = maxixs(i);
    y1 = v(rmax)-v(rmax-1);
    y2 = v(rmax+1)-v(rmax-1);
    a = (y2 - 2*y1)/2;
    b = 1-y1/a;
    xmax(i) = rmax-1+b/2;
    ymax(i) = v(rmax-1)-a*b*b/4;
  end
end
