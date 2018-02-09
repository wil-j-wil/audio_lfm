function [sinu,res,F,M,P] = sin_res(snd,fs,fft_size,win,hop)

% Calculate spectrogram and inst.frq.-gram
%fft_size = 256;
%win = 256;
%hop = 128;
[I,S]=ifgram(snd,fft_size,win,hop,fs);
% Extract the peak tracks based on the new STFT
[R,M]=extractrax(abs(S));

% Calculate the interpolated IF-gram values for exact track frequencies
F = colinterpvals(R,I);

% Interpolate the (columnwise unwrapped) STFT phases to get exact peak
% phases for every sample point (synthesis phase is negative of analysis)
P = -colinterpvals(R,unwrap(angle(S)));

% Pad each control matrix with an extra column at each end,
% to account for the N-H = H points lost fitting the first & last windows
F = [0*F(:,1),F,0*F(:,end)];
M = [0*M(:,1),M,0*M(:,end)];
P = [0*P(:,1),P,0*P(:,end)];
% (mulitplying nearest column by 0 preserves pattern of NaN values)

% only use the most prominent 3 sinusoids
%[~,mind] = max(max(M,[],2));
%F = F(1:3,:);
%M = M(1:3,:);
%P = P(1:3,:);

% Now, the phase preserving resynthesis:
sinu = synthphtrax(F,M,P,fs,win,hop)';

% Noise residual is original with harmonic reconstruction subtracted off
res = snd - sinu(1:length(snd));

end