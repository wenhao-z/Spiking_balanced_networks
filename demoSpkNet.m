% Demo of the spiking neural network 
% Wen-Hao Zhang,
% wenhao.zhang@pitt.edu
% University of Pittsburgh
% Jan 29, 2019

% Set work path
setWorkPath;

% Load the parameters
parsSpkNet;

parsNet.T = 5e3; 

parsNet.jxe = 25; % E synaptic weight
% parsNet.jxi = 40; % I synaptic weight (absolute value)
parsNet.ratiojie = 5; % The ratio between I synapse over E synapse (absolute value)
parsNet.ratio_jeestruct = 0.37;

parsNet.AmplIff = 0; %5e-4;
parsNet.IBkg_scale = [1.1; 1.05]; % The background inputs received by E and I neurons
parsNet.FanoFactorInput = 0;

% Calculate dependent parameters
parsNet = getDependentPars(parsNet);

%% Run network simulation
% Optional output arguments
outArgsOpt = struct('v', [], ...
    'bSpk', []);

tic
outSet = simSpkNet(parsNet, outArgsOpt);
toc

fireRate = accumarray(outSet.tSpk(1,:)', 1)/parsNet.T*1e3;
fprintf('Mean rate of E Neurons: %3.2fHz \n', mean(fireRate(1:parsNet.Ne)))
fprintf('Mean rate of I Neurons: %3.2fHz \n', mean(fireRate(parsNet.Ne+1:end)))

%% Plot rastergram

addpath('PlotCode')
hAxe = plotSpkRaster(outSet.tSpk, parsNet);
