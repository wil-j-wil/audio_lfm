%INT_CUBATURE Unscaled sigma points and weights for 3rd order spherical
% cubature rule
% 
% Syntax: 
%   [X,W] = int_cubature(d,dummy)
%
% In: 
%   d - State dimensionality
%
% Out:
%   X - Unscaled sigma points
%   W - Weights of sigma points
% 
% Copyright (C) 2012 Jouni Hartikainen
%
% This software is distributed under the GNU General Public 
% Licence (version 3 or later); please refer to the file 
% Licence.txt, included with the software, for details.
function [X,WM,WC] = int_cubature(d,dummy)   
    X = sqrt(d).*[eye(d) -eye(d)];
    WM = 1/(2*d).*ones(1,2*d);
    WC = WM;
end

