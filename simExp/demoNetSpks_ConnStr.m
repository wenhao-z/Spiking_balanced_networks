% Test the firing rate and its distribution with connection strength
% Author: Wen-Hao Zhang
% wenhao.zhang@pitt.edu
% University of Pittsburgh
% Feb 7, 2019

parsNet.T = 10e3;

parsNet.jxe = 15:5:40; % E synapse weight
% parsNet.jxe = [0.1, 0.2, 0.5, 1,2,5, 10:20:90]; % E synapse weight
% parsNet.jxi = 20:5:40; % I synapse weight (absolute value)
parsNet.ratiojie = 5; % The ratio between I synapse over E synapse
parsNet.ratio_jeestruct = 0.35;

parsNet.AmplIff = 0;
parsNet.IBkg_scale = [1.1; 1.05];

% Generate grid of parameters
[parGrid, dimPar] = paramGrid(parsNet);

% Calculate dependent parameters
parGrid = arrayfun(@(x) getDependentPars(x), parGrid);

%% Run network simulation
tSpkArray = cell(size(parGrid));

tStart = clock;
parfor iterPar = 1: numel(parGrid)
    fprintf('Progress: %d/%d\n', iterPar, numel(parGrid));
    outSet = simSpkNet(parGrid(iterPar));
    tSpkArray{iterPar} = outSet.tSpk;
end
tEnd = clock;

%% Perform some analysis of spikes
listAnalysis = {'PopVector', 'IntSpkIntval', 'Rate_AvgXTime', 'Rate_AvgXNeuron'};
NetStatStruct = statNetSpks([], [], listAnalysis, 'bInit', 1);
NetStatStruct = repmat(NetStatStruct, size(parGrid));

for iterPar = 1: numel(parGrid)
    NetStatStruct(iterPar) = statNetSpks(tSpkArray{iterPar}, parsNet, listAnalysis);
end

for varName = fieldnames(NetStatStruct)'
    szFields = size(NetStatStruct(1).(varName{1}));
    
    statSpk = [NetStatStruct.(varName{1})];
    statSpk = reshape(statSpk, [szFields, size(parGrid)]);
    statSpk = permute(statSpk, [ndims(szFields)+1:ndims(statSpk), 1:ndims(szFields)]); % last dim is index of group
    NetStat.(varName{1}) = squeeze(statSpk); % First several dims are the same as parGrids. The last few dims are the same as varName
end
clear varName statSpk NetStatStruct szFields iterPar


%% Plot
figure

% Mean firing rate averaged over all neurons vs. parameters
hAxe(1) = subplot(3,2,1);
hold on
% plot(dimPar.valuePar, mean(NetStat.fireRate(:,1:parsNet.Ne),2), '-b'); % E neurons
% plot(dimPar.valuePar, mean(NetStat.fireRate(:,parsNet.Ne+1:end),2), '-r'); % I neurons
errorbar(dimPar.valuePar, mean(NetStat.fireRate_neuron(:,1:parsNet.Ne),2), ...
    std(NetStat.fireRate_neuron(:,1:parsNet.Ne),0, 2), 'b'); % E neurons
errorbar(dimPar.valuePar, mean(NetStat.fireRate_neuron(:,parsNet.Ne+1:end),2), ...
    std(NetStat.fireRate_neuron(:,parsNet.Ne+1:end),0, 2), 'r'); % I neurons
ylabel('Mean firing rate (Hz)')
legend('Exc. neurons', 'Inh. neurons', 'location', 'best')
% axis square

% Firing rate distribution
hAxe(2) = subplot(3,2,2);
hold on
cSpec = cool(length(parGrid));
for iter = 1: numel(parGrid)
    rateAvg = NetStat.fireRate_neuron(iter,1:parsNet.Ne);
    [rateCount, edge] = histcounts(rateAvg, 30, 'Normalization', 'probability');
    stairs(edge(1:end-1), rateCount, 'color', cSpec(iter,:), 'linew', 1)
end
xlabel('Rate (Hz)')
ylabel('Probability')
% axis square

clear rateAvg rateCount

% CV of Inter-spike intervals
hAxe(3) = subplot(3,2,3);
hold on
CV_ISIE = reshape([NetStat.ISIRes.CV_ISIE], size(parGrid));
CV_ISII = reshape([NetStat.ISIRes.CV_ISII], size(parGrid));
plot(dimPar.valuePar, CV_ISIE, 'b')
plot(dimPar.valuePar, CV_ISII, 'r')
ylabel('CV of ISI')

% Histogram of Inter-spike intervals of E neurons
hAxe(4) = subplot(3,2,4);
hold on
for iter = 1: numel(parGrid)
    stairs(NetStat.ISIRes(iter).edgeISIE(1:end-1), NetStat.ISIRes(iter).histISIE, ...
        'color', cSpec(iter,:));
end
xlabel('Inter-spike interval (ms)')
ylabel('Probability')

% Var. of (firing rate across neurons) along time
hAxe(5) = subplot(3,2,5);
hold on
varRate_time = var(NetStat.fireRate_time, [], 3);
plot(dimPar.valuePar, varRate_time(:,1), 'b')
plot(dimPar.valuePar, varRate_time(:,2), 'r')
xlabel(dimPar.namePar)
ylabel('Var. firing rate')

linkaxes(hAxe([1,3,5]), 'x')
set(hAxe([1,3,5]), 'xscale', 'log')