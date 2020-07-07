% Calculate the statistics of spikes in network simulation
% Wen-Hao Zhang
% Math Department, University of Pittsburgh
% wenhao.zhang@pitt.edu
% Jan 30, 2019

% Run a long sequence and divide the spikes into multiple small segments

neuronEdge = [1: parsNet.Ncells, parsNet.Ncells+0.5];
tSeg = 0.5e3; % duration of each segment, unit: ms
tEdge = tSeg: tSeg: parsNet.T;

nSpk = histcounts2(outSet.tSpk(1,:)', outSet.tSpk(2,:)', neuronEdge, tEdge);


%% Factor analysis
% Get the trials with spikes
% IdxTrial = (sum(nSpk(1:parsNet.Ne,:), 1)>0);
% IdxNeuron = (sum(nSpk(1:parsNet.Ne,:),2)>0);
% 
% [LoadMat, psi] = factoran(nSpk(IdxNeuron, IdxTrial)',5);


%%
% Average firing rate over neurons
rateAvg = mean(nSpk(1:parsNet.Ne,:),2); 
fanoFactor_nSpk = var(nSpk(1:parsNet.Ne,:), 0, 2) ./ rateAvg;
rateAvg = rateAvg / tSeg * 1e3;

% Correlation coefficient of spike count between all pairs of E neurons
corrCoef = corr(nSpk(1:parsNet.Ne,:)');
corrCoef(1: size(corrCoef,1)+1:end) = 0;

%% Cross-correlation between E-E neurons (average over trials)

tBin = 5; % time bin, unit: ms
maxLag = 200; % unit: ms;
nLag = round(maxLag/tBin);
neuronEdge = 1: parsNet.Ne;
nNeurons = 50;
IdxNeuron = randperm(parsNet.Ne, nNeurons);
nTrials = round(parsNet.T/tSeg);

xCorrFun = zeros(2*nLag+1, nNeurons^2, nTrials-1);
for iter = 1: (nTrials-1)
    tEdge = 0: tBin: tSeg;
    tEdge = tEdge + iter*tSeg; % Note: throw away the 1st segment, so it is not (iter-1)*tSeg
    bSpk = histcounts2(outSet.tSpk(1,:)', outSet.tSpk(2,:)', neuronEdge, tEdge);
    [xCorrFun(:,:,iter), tLags] = xcorr(bSpk(IdxNeuron,:)', nLag);
end

tLags = tLags * tBin;
xCorrFun = permute(xCorrFun, [1,3,2]);
xCorrFun = mat2cell(xCorrFun, size(xCorrFun,1), size(xCorrFun,2), ones(1, nNeurons^2));
xCorrFun = reshape(xCorrFun, nNeurons, nNeurons);

% get all auto-correlation function
aCorrFun = xCorrFun(1:size(xCorrFun,1)+1:end);
aCorrFun = cell2mat(shiftdim(aCorrFun,-1)); % [lag, trial, cell]


xCorrFun(1:size(xCorrFun,1)+1:end) = {[]};
xCorrFun = reshape(xCorrFun, 1, []);
xCorrFun = cell2mat(shiftdim(xCorrFun,-1));
xCorrFun = reshape(xCorrFun, size(xCorrFun,1), []); % [lat, trial, cell pair]

clear bSpk tEdge

%% Plot the results
figure

% Histogram of firing rate
subplot(2,3,1)
[rateCount, edge] = histcounts(rateAvg, 30, 'Normalization', 'probability');
stairs(edge(1:end-1), rateCount)
xlabel('Rate (Hz)')
ylabel('Probability')

% Histogram of Fano factor
subplot(2,3,2)
[fanofactorCount, edge] = histcounts(fanoFactor_nSpk, 30, 'Normalization', 'probability');
stairs(edge(1:end-1), fanofactorCount)
xlabel('Fano factor')

% Mean correlation coefficient of spike count between all neuron pairs
subplot(2,3,3)
corrCoefAvg = corrCoef;
corrCoefAvg(corrCoefAvg==0) = [];
[corrCoefCount, edge] = histcounts(corrCoefAvg, 30, 'Normalization', 'probability');
stairs(edge(1:end-1), corrCoefCount)
xlabel('Corr. coef.')


% Auto-correlation function
subplot(2,3,4)
plot(tLags, mean(mean(aCorrFun,2),3))
xlabel('Lag (ms)')
ylabel('Auto correlation func.')

% Cross-correlation function
subplot(2,3,5)
plot(tLags, mean(mean(xCorrFun,2),3))
xlabel('Lag (ms)')
ylabel('Cross correlation func.')

