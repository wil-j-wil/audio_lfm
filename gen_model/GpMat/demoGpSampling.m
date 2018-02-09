%% Sample from RBF prior
gpSample('rbf', 10, [16 1], [-3 3], 1e5);

%% Sample from Matern5/2 prior
gpSample('matern52', 10, [0.8 1], [-3 3], 1e5);

%% Sample from RBF posterior
% Note: Assumes zero noise
gpPosteriorSample('rbf', 5, [16 1], [-3 3], 1e5);

%% Sample from Matern5/2 posterior with noise assumption
sig_noise = 0.25;
gpPosteriorSample('matern52', 5, [0.8 1], [-3 3], [], [], sig_noise, 1);