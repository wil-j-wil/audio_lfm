%VEC2MA Function for extracting m and A from a vector.
% 
% Syntax: 
%   [mA,m,A] = vec2mA(vmA,d)
%
% In: 
%   vmA - m and A as a vector
%
% Out:
%   mP  - Matrix containing m and A
%    m  - Mean vector
%    A  - Cholesky factor of a covariance matrix P
% 
% Copyright (C) 2012 Jouni Hartikainen
%
% This software is distributed under the GNU General Public 
% Licence (version 3 or later); please refer to the file 
% Licence.txt, included with the software, for details.
function [mA,m,A] = vec2mA(vmP,d)

    m  = vmP(1:d);
    Av = vmP(d+1:d+d*(d+1)/2);
        
    A   = zeros(d,d);
    ind = logical(tril(ones(d)));
    A(ind) = Av;
        
    mA = [m A];
end

