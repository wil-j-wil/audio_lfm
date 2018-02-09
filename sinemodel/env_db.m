function [x_rms_dBFS] = env_db(x, fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Institution: Stanford University
% Project: Thesis Research
% Author: Ryan Cassidy (05157787)
% Date: Winter 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENV_DB computes the time-domain envelope of the provided input signal, in
% dBFS.
%
% Development Status: Pre-release testing.
%
% Prior RCS Version Control Information:
% 
% $RCSfile: $
%
% $Author: $
%
% $Date: $
%
% $Locker: $
%
% $Log: $
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute RMS of x, in dBFS.
x_sq = x.^2;
tau = 20e-3; % RMS time constant seconds.
alpha_rms = exp(-1/(tau*fs));
x_ms = filter((1-alpha_rms),[1,-alpha_rms], x_sq);
x_rms = sqrt(x_ms);

x_min_dBFS = -105;

x_rms_dBFS = 20*log10(max(x_rms,10^(x_min_dBFS/20)));
