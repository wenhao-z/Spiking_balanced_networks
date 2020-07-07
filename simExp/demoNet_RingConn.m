% Demo of the spiking neural network with the ring connection structure 

% Wen-Hao Zhang,
% wenhao.zhang@pitt.edu
% University of Pittsburgh
% July 4, 2020

parsNet.T = 10e3; 

parsNet.jxe = 25; % E synaptic weight
% parsNet.jxi = 40; % I synaptic weight (absolute value)
parsNet.ratiojie = 4; % The ratio between I synapse over E synapse (absolute value)
parsNet.ratio_jeestruct = 0.3;

parsNet.AmplIff = 1e-3;
parsNet.IBkg_scale = [1; 1];
% parsNet.FanoFactorInput = 0.01;

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

%% Statistical analysis

% Projections on stimulus feature manifold
tStat = 0; % unit: ms
tBin = 5; % unit: ms
[DecodeRes, DecodePars] = decoder_PopVector(outSet.tSpk, parsNet, 'tStat', tStat, ...
    'tBin', tBin, 'bOutBumpPos', 1);

% Cross correlation function
BumpPos = DecodeRes.BumpPos;
BumpPos(isnan(BumpPos)) = 0;
[CCFunc, tLag] = xcorr(BumpPos, 50);
CCFunc = CCFunc((end+1)/2:end)./max(CCFunc);
tLag = tLag((end+1)/2:end) * tBin;

clear BumpPos
%% Plot 

addpath('PlotCode')
% rastergram
hAxe = plotSpkRaster(outSet.tSpk, parsNet);

%%
figure;
hAxe2(1) = subplot(2,4,1:3);
plot(tBin: tBin: parsNet.T, DecodeRes.BumpPos)
xlabel('Time (ms)')
ylabel('Stimulus sample')

hAxe2(2) = subplot(2,4,4);
hold on
[histBumpPos, edgeSample] = histcounts(DecodeRes.BumpPos, 100);
stairs(histBumpPos/sum(histBumpPos), (edgeSample(1:end-1)+edgeSample(2:end))/2);
plot(normpdf((edgeSample(1:end-1)+edgeSample(2:end))/2, DecodeRes.meanBumpPos, sqrt(DecodeRes.varBumpPos))*mean(diff(edgeSample)), ...
    (edgeSample(1:end-1)+edgeSample(2:end))/2)
xlabel('Dist.of samples')

linkaxes(hAxe2, 'y')

% Cross correlation function
subplot(2,2,3)
plot(tLag, CCFunc)
ylabel('Cross-correlation')
xlabel('Time (ms)')