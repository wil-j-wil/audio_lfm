%% Generates novel instances of the natural sounds recordings

% Will Wilkinson Jan-2018

% include_ind selects just the sounds that we used for the experiment
include_ind = [1:12 14 15 18 19 21:24 26];
sig_duration = 4; % in seconds
suffix = {'_1','_2'};
for i=include_ind %[3 8 22 26]
    [cascadeModLFM, genLFM] = genModelLFM(i,sig_duration,suffix); % Generate with LFM gen model
    [cascadeModNMF, genNMF] = genModelNMF(i,sig_duration,suffix); % Generate with NMF gen model
    [cascadeModtNMF, gentNMF] = genModelTNMF(i,sig_duration,suffix); % Generate with tNMF gen model
    genAnchor(i,sig_duration); % Generate anchor
end