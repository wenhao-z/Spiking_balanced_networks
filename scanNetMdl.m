% Scan the parameters for the spiking neural network
% Wen-Hao Zhang,
% wenhao.zhang@pitt.edu
% University of Pittsburgh
% Feb. 7, 2019

% Set work path
setWorkPath;

% Load the parameters
parsSpkNet;

parsNet.T = 20e3;
parsNet.jxe = 15:5:50; % E synapse weight
parsNet.jxi = 125; % E synapse weight
parsNet = rmfield(parsNet, 'jxi');
% parsNet.ratiojie = 4:1:7; % The ratio between I synapse over E synapse
parsNet.ratio_jeestruct = 0.1:0.1:0.4;

parsNet.AmplIff = [0.1, 0.2, 0.5, 1,2,5]*1e-2;
parsNet.FanoFactorInput = [0, 0.01, 0.05, 0.1];

% Generate grid of parameters
[parGrid, dimPar] = paramGrid(parsNet);

% Calculate dependent parameters
parGrid = arrayfun(@(x) getDependentPars(x), parGrid);

%% Run network simulation
tSpkArray = cell(size(parGrid));

parpool(12);
tStart = clock;
parfor iterPar = 1: numel(parGrid)
    fprintf('Progress: %d/%d\n', iterPar, numel(parGrid));    
    outSet = simSpkNet(parGrid(iterPar));
    tSpkArray{iterPar} = outSet.tSpk;
end
tEnd = clock;

% Conver the format of tSpk in order to save space
% tSpkArray = cellfun(@(x) tSpkFormatConverter(x, parsNet.Ncells), tSpkArray, 'uniformout', 0);

%% Save
savePath = fullfile(Path_RootDir, 'Data');
mkdir(savePath);

str = datestr(now, 'yymmddHHMM');
fileName = ['scanSpkNetPars_', str(1:6), ...
    '_', str(7:end) '.mat'];

% save(fullfile(savePath, fileName))
save(fullfile(savePath, fileName), '-v7.3')
