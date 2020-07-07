% Demo of the spiking neural network through plotting the trajectory of
% membrane potential and rastergram.

% Wen-Hao Zhang,
% wenhao.zhang@pitt.edu
% University of Pittsburgh
% Feb 9, 2019

parsNet.T = 2e3; 

parsNet.jxe = 25; % E synaptic weight
% parsNet.jxi = 40; % I synaptic weight (absolute value)
parsNet.ratiojie = 4; % The ratio between I synapse over E synapse (absolute value)
parsNet.ratio_jeestruct = 0.3;

% parsNet.IBkg_scale = 0.8*[1; 1];
% parsNet.FanoFactorInput = 0.01;

parsNet.AmplIff = 0;
% Calculate dependent parameters
parsNet = getDependentPars(parsNet);

%% Run network simulation
% Optional output arguments
outArgsOpt = struct('v', [], ...
    'bSpk', []);

tic
outSet = simSpkNet_vtrace(parsNet, outArgsOpt);
toc

fireRate = accumarray(outSet.tSpk(1,:)', 1)/parsNet.T*1e3;
fprintf('Mean rate of E Neurons: %3.2fHz \n', mean(fireRate(1:parsNet.Ne)))
fprintf('Mean rate of I Neurons: %3.2fHz \n', mean(fireRate(parsNet.Ne+1:end)))

%% Plot 

addpath('PlotCode')
% rastergram
hAxe = plotSpkRaster(outSet.tSpk, parsNet);

numAxes = length(hAxe);
tStat = 100; %unit:ms
nBins = 300;

IdxNeuronE = randperm(parsNet.Ne, 1);
IdxNeuronI = randperm(parsNet.Ni, 1) + parsNet.Ne;
tSpkE = outSet.tSpk(2, outSet.tSpk(1,:)==IdxNeuronE);
tSpkI = outSet.tSpk(2, outSet.tSpk(1,:)==IdxNeuronI);

% Mark the selected neurons in rastergram
axes(hAxe(1)); hold on
scatter(tSpkE, parsNet.PrefStim(IdxNeuronE)*ones(size(tSpkE)), 'om')
axes(hAxe(2)); hold on
scatter(tSpkI, IdxNeuronI*ones(size(tSpkI))-parsNet.Ne, 'og')

% Plot membrane potenatial trace
figure; 

hAxe(numAxes+1) = subplot(2,5,(1:4));
hold on
tPlot = (0:parsNet.T/parsNet.dt) * parsNet.dt;
% Membrane potential trace of two example neurons
plot(tPlot, outSet.v(IdxNeuronE,:), 'b');
plot(tPlot, outSet.v(IdxNeuronI,:), 'r');
% Corresponding spike timing
plot(repmat(tSpkE, 2,1), repmat([1.1; 1.2], 1, length(tSpkE)), 'b');
plot(repmat(tSpkI, 2,1), repmat([1.1; 1.2], 1, length(tSpkI)), 'r');
title(sprintf('Neuron index (E:%d, I:%d).', IdxNeuronE, IdxNeuronI))

hAxe(numAxes+2) = subplot(2,5,(1:4)+5);
hold on
plot(tPlot, mean(outSet.v(1:parsNet.Ne,:),1), 'b');
plot(tPlot, mean(outSet.v(parsNet.Ne+1:end,:)), 'r');
title(sprintf('Average over neurons'))
xlabel('Time (ms)')
ylabel('Membrane potential (mV)')

hAxe(numAxes+3) = subplot(2,5,5); 
hold on
% Distribution of membrane potential in equilibrium state
[histVE, vEdgeE] = histcounts(outSet.v(tStat/parsNet.dt+1, 1:parsNet.Ne), nBins, 'normalization', 'probability');
[histVI, vEdgeI] = histcounts(outSet.v(tStat/parsNet.dt+1, parsNet.Ne+1:end), nBins, 'normalization', 'probability');
stairs(histVE, (vEdgeE(1:end-1)+vEdgeE(2:end))/2, 'b');
stairs(histVI, (vEdgeI(1:end-1)+vEdgeI(2:end))/2, 'r');

linkaxes(hAxe(1:5), 'x')
linkaxes(hAxe([4,6]), 'y')

