% MYODE Wrapper function used with ODE solvers in moment integration
% 
% Syntax: 
%   df = myode(t,x,func,param)
%
% In: 
%    t     - Time instance
%    x     - State vector
%    func  - Moment derivative function with interface ()
%    param - Parameters of moment derivative function
%
% Out:
%       df - Derivatives of moments
% 
% Copyright (C) 2012 Jouni Hartikainen
%
% This software is distributed under the GNU General Public 
% Licence (version 3 or later); please refer to the file 
% Licence.txt, included with the software, for details.
function [df] = myode(t,x,func,param)
    param{end} = t;
    df = feval(func,x,param);
end

