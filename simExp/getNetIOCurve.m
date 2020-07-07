% Get the input-firing rate curve of network model
% Author: Wen-Hao Zhang
% wenhao.zhang@pitt.edu
% University of Pittsburgh
% Feb 9, 2019

parsNet.T = 4e3;

parsNet.jxe = 25; % E synapse weight
parsNet.ratiojie = 5; % The ratio between I synapse over E synapse
parsNet.ratio_jeestruct = 0;

parsNet.AmplIff = 0;
% parsNet.IBkg_scale = [0.9, 0.9];

parsNet.IBkg_scale = [1; 1];
parsNet.IBkg_scale = parsNet.IBkg_scale * (0.4:0.1:1); %(0: 0.3: 1.5);
parsNet.IBkg_scale = parsNet.IBkg_scale + 0*ones(2,1);
% (0.4:0.1:1), (1: 0.2: 2.4)

% gain = 0.5: 0.1: 1.1;
% gain = [gain; ones(size(gain))];
% parsNet.IBkg = bsxfun(@times, parsNet.IBkg, gain);

parsNet.FanoFactorInput = [0.1, 0.5]; %[0.01, 0.05];

% Generate grid of parametersx
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
listAnalysis = {'Rate_AvgXTime', 'fanoFactor_nSpk'}; 
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


%% Theoretical prediction of mean firing rate

addpath(fullfile(Path_RootDir, 'TheoreticalAnalysis'));
rateTheory = arrayfun(@(netpars) getFiringRateTheory(netpars), shiftdim(parGrid, -1), 'uniformout', 0);
rateTheory = cell2mat(rateTheory); % [dim of results, size(parGrid)]


%% Plot code

dim_avg = length(dimPar)+1;
IdxDim_IBkg = cellfun(@(x) strcmp(x, 'IBkg_scale'), {dimPar.namePar});
IdxDim_IBkg = find(IdxDim_IBkg==1);
fireRate_neuron = permute(NetStat.fireRate_neuron, [dim_avg, 1:dim_avg-1]);
fanoFactor_nSpk = permute(NetStat.fanoFactor_nSpk, [dim_avg, 1:dim_avg-1]);

cSpec = cool(length(parsNet.FanoFactorInput));
% cSpec = cool(length(dimPar));
% cSpecI = hot(length(dimPar));

figure;
% Plot the firing rate vs. input current
hAxe(1) = subplot(2,2,1);
hold on
for iterFF = 1: length(parsNet.FanoFactorInput)
    plot(dimPar(IdxDim_IBkg).valuePar(1,:), ...
        squeeze(mean(fireRate_neuron(1:parsNet.Ne,iterFF,:),1)), ...
        'o', 'color', cSpec(iterFF,:)); % E neurons
    plot(dimPar(IdxDim_IBkg).valuePar(1,:), squeeze(rateTheory(1,iterFF,:)), ...
        'color', cSpec(iterFF,:)); % E neurons, theoretical prediction
end
% axis square
title('Exc neurons')
xlabel(dimPar(IdxDim_IBkg).namePar)
ylabel('Firing rate (Hz)')

hAxe(2) = subplot(2,2,2);
hold on
for iterFF = 1: length(parsNet.FanoFactorInput)
    plot(dimPar(IdxDim_IBkg).valuePar(1,:), ...
        squeeze(mean(fireRate_neuron(parsNet.Ne+1:end,iterFF,:),1)), ...
        's', 'color', cSpec(iterFF,:)); % I neurons
    plot(dimPar(IdxDim_IBkg).valuePar(1,:), squeeze(rateTheory(2,iterFF,:)), ...
        'color', cSpec(iterFF,:)); % I neurons, theoretical prediction
end
% axis square
title('Inh neurons')
xlabel(dimPar(IdxDim_IBkg).namePar)
ylabel('Firing rate (Hz)')

% Plot the Fano factor vs. input current
hAxe(3) = subplot(2,1,2);
hold on
for iterFF = 1: length(parsNet.FanoFactorInput)
    plot(dimPar(IdxDim_IBkg).valuePar(1,:), ...
        squeeze(mean(fanoFactor_nSpk(1:parsNet.Ne,iterFF,:),1)), ...
        'color', cSpec(iterFF,:)); % E neurons
    plot(dimPar(IdxDim_IBkg).valuePar(1,:), ...
        squeeze(mean(fanoFactor_nSpk(parsNet.Ne+1:end,iterFF,:),1)), ...
        '--', 'color', cSpec(iterFF,:)); % I neurons
end
% axis square
xlabel(dimPar(IdxDim_IBkg).namePar)
ylabel('Fano factor')

linkaxes(hAxe, 'x')