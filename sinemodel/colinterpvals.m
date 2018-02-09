function V = colinterpvals(T,M)
% V = colinterpvals(T,M)  Interpolate columns from real-valued bin tracks
%      Each row of T defines a 'track' as a continuously-valued 'bin index'
%      (where zero indicates the 1st bin i.e. NOT matlab-style indices)
%      for the corresponding column of matrix M; the conformal value of 
%      V is the linear interpolation of M at that bin.
% 1998may03 dpwe@icsi.berkeley.edu for AHI tracking 
% $Header: $

[nr, nc] = size(M);

% Preallocate results array
V = zeros(size(T));

[nt,ntc] = size(T);

for c = 1:nc
  mc = M([[1:nr] nr],c);
  tc = T(:,c)+1;  % add 1 to make it valid matlab indexing
  bin = floor(tc);
  binOK = (bin >= 1) & (bin <= nr);
  bin(find(~binOK)) = ones(1,sum(~binOK));
  base = mc(bin);
  delta = mc(bin+1) - base;
  % adding the tc restores the NaNs
  vals = base+(tc-bin).*delta;
  V(:,c) = vals;
end
