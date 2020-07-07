function [hFig, hAxe] = plotDecodeResFunc(plotDatStruct, IdxPars, dimSpecPlot, hFig)
% Plot the statistics of bump position with network parameters
% Wen-Hao Zhang,
% wenhao.zhang@pitt.edu
% University of Pittsburgh
% Jan 29, 2019

% plotDataStruct: .....
% IdxPars: the index of the example parameter in parGrid
%% Unfold the fields from struct
for varName = fieldnames(plotDatStruct)'
    eval([varName{1} '=plotDatStruct.' varName{1} ';']);
end

nPars = length(dimPar);

%% Generate Fig. and Axes handle
scrsz = get(0,'ScreenSize');
% close all;
if isempty('hFig')
    hFig = figure('Position',[scrsz(4)/3 scrsz(4)/4 scrsz(3)*.6 scrsz(4)*.6]);
else
    set(hFig, 'Position',[scrsz(4)/3 scrsz(4)/4 scrsz(3)*.6 scrsz(4)*.6]);
    clf;
end

for iter = 1: 2*nPars
    hAxe(iter) = subplot(2, nPars,iter);
    axis square;
    hold on
end


%% Mean firing arte
Rate = permute(NetStat.fireRate_neuron, [length(dimPar)+1, 1:length(dimPar)]);
szRate = size(Rate);
Rate = reshape(Rate, szRate(1), []);

NetStat.meanRateE = squeeze(mean(Rate(1:parsNet.Ne,:),1));
NetStat.meanRateI = squeeze(mean(Rate(parsNet.Ne+1:end,:),1));

NetStat.meanRateE = reshape(NetStat.meanRateE, szRate(2:end));
NetStat.meanRateI = reshape(NetStat.meanRateI, szRate(2:end));
clear Rate szRate

%% 
% Get the index of example parameters 
IdxParsExamp = zeros(1, length(dimPar));
for varName = fieldnames(IdxPars)'
    idx = cellfun(@(x) strcmp(x, varName{1}), {dimPar.namePar});
    IdxParsExamp(idx) = IdxPars.(varName{1});   
end

for iterPar = 1: nPars
    % Loop for the parameter shown on the x-axis   
    
    for varName = {'concBumpPos', 'meanRateE', 'meanRateI'}
        plotDat = NetStat.(varName{1});
        szDat = size(plotDat);
        cumprodGrid = cumprod(szDat);
        
        % Calculate the index of dataum to be plotted out in plotDat 
        nValues = length(dimPar(iterPar).valuePar);
        subs = repmat(IdxParsExamp, [nValues, ones(1, nPars-1)]);
        subs(:,iterPar) = 1: nValues;
        
        IdxDat = (subs(:, 2:end)-1) .* cumprodGrid(1:end-1);
        IdxDat = sum(IdxDat,2) + subs(:,1);
        
        switch varName{1}
            case 'concBumpPos'
                plot(hAxe(iterPar+nPars), dimPar(iterPar).valuePar, plotDat(IdxDat), '-o');
            otherwise
                plot(hAxe(iterPar), dimPar(iterPar).valuePar, plotDat(IdxDat), '-o');
        end
    end   
    xlabel(hAxe(iterPar+nPars), dimPar(iterPar).namePar)
end

ylabel(hAxe(1+nPars), 'Conc. of bump pos.')
ylabel(hAxe(1), 'Mean firing rate (Hz)')
legend(hAxe(1), 'Exc. neurons', 'Inh. neurons', 'location', 'best')

strTitle = cell(1, length(dimPar));
for iterPar = 1: length(dimPar)
    namePar = dimPar(iterPar).namePar;
    strTitle(iterPar) = { [namePar , '=', ...
        num2str(dimPar(iterPar).valuePar(IdxPars.(namePar)))]};
end
title(hAxe(2), strTitle)

