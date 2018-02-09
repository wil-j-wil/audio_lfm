% MA2VEC Converts matrix containing mean and chol of covariance to a vector
% 
% Syntax: 
%   [vmA,m,A] = mA2vec(mA,d)
%
% In: 
%   mA - Matrix containing mean and chol of a covariance matrix
%    d - Length of mean vector
%
% Out:
%   vmA - mA in vector form
%     m - Mean vector
%     A - Cholesky of a covariance matrix
% 
% Copyright (C) 2012 Jouni Hartikainen
%
% This software is distributed under the GNU General Public 
% Licence (version 3 or later); please refer to the file 
% Licence.txt, included with the software, for details.


function [vmA,m,A] = mA2vec(mA,d)
        
    m = mA(1:d);
    A = mA(d+1:d+d^2);
    A = A(logical(tril(ones(d))));
    vmA = [m(:); A(:)];

end

