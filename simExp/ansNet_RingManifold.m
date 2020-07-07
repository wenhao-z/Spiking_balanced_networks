% Demo of the spiking neural network
% Wen-Hao Zhang,
% wenhao.zhang@pitt.edu
% University of Pittsburgh
% Jan 29, 2019


parsNet.T = 10e3;

% Fix E synapse, and increase I synapse
% parsNet.jxe = 25; % E synapse weight
% parsNet.ratiojie = 5:8; % The ratio between I synapse over E synapse

parsNet.jxe = 15:5:40; % E synapse weight
parsNet.jxi = 125; % E synapse weight
parsNet = rmfield(parsNet, 'jxi');
% parsNet.ratiojie = 5; % The ratio between I synapse over E synapse
parsNet.ratio_jeestruct = 0.3;

parsNet.AmplIff = 1e-2/2; % 1e-2;
parsNet.IBkg_scale = [1.1; 1.05]; % The background inputs received by E and I neurons
parsNet.FanoFactorInput = 0;

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

%% Calculate statistics of network spikes
tBin = 100; % May need to use a long time bin to cancel the influence from 
%            oscillation in estimating the variance of bump position

listAnalysis = {'PopVector', 'Rate_AvgXTime'};
NetStatStruct = statNetSpks([], [], listAnalysis, 'bInit', 1);
NetStatStruct = repmat(NetStatStruct, size(parGrid));

for iterPar = 1: numel(parGrid)
    NetStatStruct(iterPar) = statNetSpks(tSpkArray{iterPar}, parsNet, listAnalysis, 'tBin', tBin);
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
% addpath('PlotCode')
% hAxe = plotSpkRaster(tSpkArray{4}, parsNet);

figure
% Firing rate of all neurons
subplot(2,3,1)
hold on
% cSpec = cool(length(parGrid));
cSpec = lines(length(parGrid));
for iter = 1: length(dimPar.valuePar)
    plot(parsNet.PrefStim, smooth(NetStat.fireRate_neuron(iter,1:parsNet.Ne)', 5), ...
        'color', cSpec(iter,:))
end
% plot(parsNet.PrefStim, NetStat.fireRate_neuron(:,1:parsNet.Ne)')
xlim(parsNet.PrefStim([1,end]))
set(gca, 'xtick', parsNet.PrefStim([1,end/2,end]), 'xticklabel', parsNet.Width*[-1,0,1])
xlabel('Neuron index')
ylabel('Firing rate (Hz)')
axis square

% Mean firing rate averaged over all neurons vs. parameters
subplot(2,3,4)
hold on
plot(dimPar.valuePar, mean(NetStat.fireRate_neuron(:,1:parsNet.Ne),2), '-o'); % E neurons
plot(dimPar.valuePar, mean(NetStat.fireRate_neuron(:,parsNet.Ne+1:end),2), '-o'); % I neurons
xlabel(dimPar.namePar)
ylabel('Mean firing rate (Hz)')
legend('Exc. neurons', 'Inh. neurons', 'location', 'best')
axis tight; axis square

% Concentration of bump position vs. parameters
subplot(2,3,5)
hold on
plot(dimPar.valuePar, NetStat.concBumpPos, '-o');
xlabel(dimPar.namePar)
ylabel('Conc. of bump position')
axis tight; axis square

% Concentration of bump position vs. firing rate
subplot(2,3,6)
plot(mean(NetStat.fireRate_neuron(:,1:parsNet.Ne),2), NetStat.concBumpPos, '-o')
xlabel('Mean firing rate (Hz)')
ylabel('Conc. of bump position')
axis tight; axis square

IdxPars.AmplIff = 4; %4;
IdxPars.FanoFactorInput = 2;
IdxPars.jxe = 3; %3;
IdxPars.ratio_jeestruct = 3;

str = {['\alpha_{ff}=', num2str(parsNet.AmplIff)], ...
    ['F=' num2str(parsNet.FanoFactorInput)], ...
['\beta_{jeestruct}=' num2str(parsNet.ratio_jeestruct)]};

subplot(2,3,2)
text(0.5, 0.5, str)
axis off