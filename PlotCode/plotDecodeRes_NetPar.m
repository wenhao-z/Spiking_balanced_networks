% Plot the network estimation results with one parameter
% Wen-Hao Zhang, April-8, 2016

%% Load the data
bSavePlot = 0;

setWorkPath;
datPath = fullfile(Path_RootDir, 'Data');
addpath(fullfile(Path_RootDir, 'PlotCode'));

% fileName = 'scanSpkNetPars_190401_1740.mat';
fileName = 'BumpStatNetPars_190402_1911.mat';

plotDatStruct = load(fullfile(datPath, fileName), ...
    'NetStat', 'dimPar', 'parsNet');

%% Parameters of plotted data
dimSpecPlot.xValDimName = 'jxe';
dimSpecPlot.colorDimName = 'FanoFactorInput';
% dimSpecPlot.markerDimName = 'IdxNeuronGroup';

% The example parameters
IdxPars.AmplIff = 3; %4;
IdxPars.FanoFactorInput = 1;
IdxPars.jxe = 3; %3;
IdxPars.ratio_jeestruct = 3;


[hFig, hAxe] = plotDecodeResFunc(plotDatStruct, IdxPars, dimSpecPlot, figure(2));

%%
if bSavePlot
    cd([datPath, '/figure']);
    set(gcf, 'PaperOrientation', 'landscape')
    set(gcf, 'PaperPosition', [0.63, 0.63, 28.41, 19.72])
    saveas(gcf, 'CANN_MeanvsNoise.eps', 'psc')
end