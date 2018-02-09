%VEC2MAC Function for extracting m, A and C from a vector.
% 
% Syntax: 
%   [mAC,m,A,C] = vec2mAC(vmAC,d)
%
% In: 
%   vmAC - m, A and C as a vector
%
% Out:
%   mPC - Matrix containing m, A and C
%     m - Mean vector
%     A - Cholesky factor of a covariance matrix P
%     C - Cross covariance matrix used in smoothing
% 
% Copyright (C) 2012 Jouni Hartikainen
%
% This software is distributed under the GNU General Public 
% Licence (version 3 or later); please refer to the file 
% Licence.txt, included with the software, for details.
function [mAC,m,A,C] = vec2mAC(vmAC,d)

    m  = vmAC(1:d);
    Av = vmAC(d+1:d+d*(d+1)/2);
    C  = vmAC(d+1+d*(d+1)/2:d+d*(d+1)/2+d^2);
    
    m = reshape(m,d,1);
    
    A = zeros(d,d);
    ind = logical(tril(ones(d)));
    A(ind) = Av;
    C = reshape(C,d,d);
    
    mAC = [m A C];
end

