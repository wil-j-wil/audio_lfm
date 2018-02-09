function [X,ph,pp] = synthphtrax(F, M, P, SR, WIN, HOP, DUR)
% X = synthphtrax(F, M, P, SR, WIN, HOP, DUR) Reconstruct from tracks incl phase.
%	Each row of F, M and P contains a series of frequency, magnitude 
%       and phase samples for a particular track.  These will be reconstructed
%	and summed into the output sound X which will run at sample rate SR, 
%	although the columns in F and M are derived from an STFT with window 
%   WIN and advance HOP.
%   If DUR is nonzero, X will be padded or truncated 
% 	to correspond to just this much time.
% dpwe@icsi.berkeley.edu 1994aug20, 1996aug22

if(nargin<7)
  DUR = 0;
end

rows = size(F,1);
cols = size(F,2);

cols = size(F,2);

opsamps = round(DUR*SR);
if(DUR == 0)
  opsamps = WIN+((cols-1)*HOP);
end

X = zeros(1, opsamps);

% Time step
dt = HOP/SR;

% quadratic or cubic interpolation
iorder = 3;

for row = 1:rows

  mm = M(row,:);
  ff = F(row,:);
  pp = P(row,:);
  % Where mm = 0, ff is undefined.  But interp will care, so set them.

  % First, set all nan values of mm to zero
  mm(find(isnan(mm))) = zeros(1, sum(isnan(mm)));
  ff(find(isnan(ff))) = zeros(1, sum(isnan(ff)));
  pp(find(isnan(pp))) = zeros(1, sum(isnan(pp)));

  % For speed, chop off regions of initial and final zero magnitude - 
  % but want to include one zero from each end if they are there 
  nzv = find(ff);
  firstcol = min(nzv);
  lastcol = max(nzv);
  zz = [max(1, firstcol-1):min(cols,lastcol+1)];
  mm = mm(zz);
  ff = ff(zz);
  pp = pp(zz);
  nzcols = length(zz);

  % Find onsets - points where mm goes from zero (or NaN) to nzero
  mz = (ff==0);
  mask = mz & (0==[mz(2:nzcols),1]);
  % Set frequencies at bins before onsets to match bin at onset
  ff = ff.*(1-mask) + mask.*[ff(2:nzcols),0];
  % Wind phase back to exactly match that freq
  pp = pp.*(1-mask) + mask.*[pp(2:nzcols)-ff(2:nzcols)*dt*2*pi,0];
  % Do offsets too
  mask = mz & (0==[1,mz(1:(nzcols-1))]);
  ff = ff.*(1-mask) + mask.*[0,ff(1:(nzcols-1))];
  pp = pp.*(1-mask) + mask.*[0,pp(1:(nzcols-1))+ff(1:(nzcols-1))*dt*2*pi];
  
  % Ok. Can interpolate now
  % This is actually the slow part

  ph = zeros(1,1+(nzcols-1)*HOP);

  for col = 1:(nzcols-1)
    % time index
    tt = [0:(HOP-1)]/HOP*dt;
    dp = pp(col+1) - pp(col);
    ww = 2*pi*ff(col+[0 1]);
    % How many cycles at given frq?
    ncycu = ((ww(1) + ww(2)) * dt/2 - dp)/(2*pi);
    ncyc = round( ncycu );
    dp = dp + 2*pi*ncyc;
    if iorder == 2
      % Re-estimate final frq to match final phase
      ww(2) = 2/dt * dp - ww(1);
      % .. and recover quadratic phase coefficient
      a = (ww(2) - ww(1))/(2*dt);
      % Hence phase is...
      ph((col-1)*HOP + [1:HOP]) = pp(col) + ww(1)*tt + a*(tt.^2);
    else
      % iorder = 3
      a2 = 3/(dt*dt)*(dp - dt/3*(2*ww(1)+ww(2)));
      a3 = 1/(3*dt*dt)*(ww(2)-ww(1)-2*a2*dt);
      ph((col-1)*HOP + [1:HOP]) = pp(col) + ww(1)*tt + a2*(tt.^2) + a3*(tt.^3);
    end    
  end

  % Last value in phase function directly from input
  ph((nzcols-1)*HOP + 1) = pp(nzcols);

  % Simple linear interpolation of magnitudes
  mmi = slinterp(mm, HOP);
  
  % run the oscillator and apply the magnitude envelope
  xx = mmi.*cos(ph);  

  % overlap-add it in to the correct place in the array
  base = 1+WIN/2-HOP+HOP*(zz(1)-1);
  sizex = prod(size(xx))-1;
  xx = xx(1:sizex);
  ww = (base-1)+[1:sizex];
  X(ww) = X(ww) + xx;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Y = slinterp(X,F)
% Y = slinterp(X,F)  Simple linear-interpolate X by a factor F
%        Y will have ((size(X)-1)*F)+1 points i.e. no extrapolation
% dpwe@icsi.berkeley.edu  fast, narrow version for SWS

% Do it by rows

sx = prod(size(X));

% Ravel X to a row
X = X(1:sx);
X1 = [X(2:sx),0];

XX = zeros(F, sx);

for i=0:(F-1)
  XX((i+1),:) = ((F-i)/F)*X + (i/F)*X1;
end

% Ravel columns of X for output, discard extrapolation at end
Y = XX(1:((sx-1)*F+1));
