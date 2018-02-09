%KF_UPDATE  Square Root Kalman Filter update step
%
% Syntax:
%   [X,S,K,IM,P,LLH] = SRKF_UPDATE(X,S,Y,H,sR)
%
% In:
%   X - Nx1 mean state estimate after prediction step
%   S - NxN square root of a state covariance after prediction step
%   Y - Dx1 measurement vector.
%   H - Measurement matrix.
%   sR - Measurement noise covariance.
%
% Out:
%   X  - Updated state mean
%   P  - Update Covariance matrix of state
%   S  - Updated square root of a state covariance matrix
%   K  - Computed Kalman gain
%   LLH - Predictive log probability (likelihood) of measurement
%   
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function [X,P,S,K,LLH] = srkf_update(X,S,Y,H,sR)

  %
  % Check which arguments are there
  %
  if nargin < 5
    error('Too few arguments');
  end

  %
  % update step
  %
  IM = H*X;
  
  D = H*S;
    
  [tmp,Tr] = qr([D sR; S zeros(size(S,1),size(sR,2))]',0);
  Tr = Tr';
  T11 = Tr(1:size(Y,1),1:size(Y,1));
  T21 = Tr(size(Y,1)+1:end,1:size(Y,1));
  T22 = Tr(size(Y,1)+1:end,size(Y,1)+1:end);
  K = T21/T11;
  
  res = Y - IM;
  X = X + K*res;
  S = T22;
  P = S*S';

  if nargout > 4
    dy = size(Y,1);
    b=T11\res;
    E=-.5*b'*b;
    E=E-sum(log(diag(T11)))-.5*dy*log(2*pi);
    LLH = real(E);
  end
  

  

  
