% MA2VEC Converts matrix of mean, chol of cov. and cross cov to a vector
% 
% Syntax: 
%   [vmA,m,A,C] = mAC2vec(mA,d)
%
% In: 
%   mA - Matrix containing mean, chol of a cov matrix and cross cov
%    d - Length of mean
%
% Out:
%   vmAC - mA in vector form
%      m - Mean vector
%      A - Cholesky of a covariance matrix
%      C - Cross covariance matrix
% 
% Copyright (C) 2012 Jouni Hartikainen
%
% This software is distributed under the GNU General Public 
% Licence (version 3 or later); please refer to the file 
% Licence.txt, included with the software, for details.
function [vmAC,m,A,C] = mAC2vec(mAC,d)
        
    m = mAC(1:d);
    A = mAC(d+1:d+d^2);
    A = A(logical(tril(ones(d))));
    C = mAC(d+1+d^2:d+2*d^2);
    vmAC = [m(:); A(:); C(:)];

end

