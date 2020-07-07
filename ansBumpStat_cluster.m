function ansBumpStat_cluster(fileName)

% Load tSpkArray under different combinations of network parameters
% Analyze the statistics of bump activity.

% Wen-Hao ZHang
% wenhao.zhang@pitt.edu
% University of Pittsburgh
% Apr 3, 2019


% Load the data
setWorkPath;
addpath(fullfile(Path_RootDir, 'SpikeAnalysisKit'));

datPath = fullfile(Path_RootDir, 'Data');

fprintf(['Loading ' fileName, '.\n']);
load(fullfile(datPath, fileName));
fprintf('Loading complete.\n');

%% Calculate statistics of network spikes

fprintf('Analyze the bump activity.\n');

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

%% Save
fprintf('Save analyzed results.\n');

Idx = strfind(fileName, '_');
fileName = ['BumpStatNetPars', fileName(Idx:end)];

clear tSpkArray
save(fullfile(savePath, fileName), '-v7.3')

end