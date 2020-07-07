% Demo of the chaotic state in balanced spiking neural network
% through comparing the network activity with/without deleting only ONE SPIKE
% Wen-Hao Zhang,
% wenhao.zhang@pitt.edu
% University of Pittsburgh
% Jan 30, 2019

parsNet.T = 2e3;

parsNet.jxe = 25;
parsNet.ratio_jeestruct = 0;

% Calculate dependent parameters
parsNet = getDependentPars(parsNet);
%% Run network simulation
% Optional output arguments
outArgsOpt = struct('v', [], ...
    'bSpk', []);

outSet = simSpkNet(parsNet, outArgsOpt);

% Randomly delete one spike and continue running network dynamics
bSpk = outSet.bSpk;
IdxSpk = find(bSpk);
IdxSpkDel = randperm(length(IdxSpk), 1);
bSpk(IdxSpk(IdxSpkDel)) = 0;

% Original state
outSet0 = simSpkNet(parsNet, outArgsOpt, 'v', outSet.v, 'bSpk', outSet.bSpk);
% Delete one spike
outSet1 = simSpkNet(parsNet, outArgsOpt, 'v', outSet.v, 'bSpk', bSpk);

tSpk0 = outSet0.tSpk;
tSpk0(2,:) = tSpk0(2,:) + parsNet.T;
tSpk0 = [outSet.tSpk, tSpk0];

tSpk1 = outSet1.tSpk;
tSpk1(2,:) = tSpk1(2,:) + parsNet.T;

%% Plot
figure;
hold on

% ------------------------------------------------------
% Color for the original sequence and altered sequence
cSpecE = 'bm';
cSpecI = 'rg';

hAxe(1) = subplot(6,1, 1:4);
hold on
hAxe(2) = subplot(6,1, 5);
hold on
hAxe(3) = subplot(6,1,6);
hold on

namePlotVar = {'tSpk0', 'tSpk1'};
for iter = 1: length(namePlotVar)
    tSpk = eval(namePlotVar{iter});
    
    % Rastergram of E neurons
    IdxSpkE = (tSpk(1,:) <= parsNet.Ne);
    scatter(hAxe(1), tSpk(2,IdxSpkE), parsNet.PrefStim(tSpk(1,IdxSpkE)), 1, cSpecE(iter))
    
    % Rastergram of I neurons
    IdxSpkI = (tSpk(1,:) > parsNet.Ne);
    scatter(hAxe(2), tSpk(2,IdxSpkI), tSpk(1,IdxSpkI)-parsNet.Ne, 1, cSpecI(iter))
    
    % ------------------------------------------------------
    % Mean firing rate of E and I neurons across time   
    if ~exist('tBin', 'var')
        tBin = 5; % unit: ms
    end
    tEdge = [0: tBin : 2*parsNet.T, 2*parsNet.T + tBin/2];
    
    neuronEdge = [1, parsNet.Ne + 0.5];
    bSpkE = histcounts2(tSpk(1,:)', tSpk(2,:)', neuronEdge, tEdge);
    bSpkE = bSpkE/ parsNet.Ne/ tBin *1e3;
    
    neuronEdge = [parsNet.Ne+1, parsNet.Ncells+0.5];
    bSpkI = histcounts2(tSpk(1,:)', tSpk(2,:)', neuronEdge, tEdge);
    bSpkI = bSpkI/ parsNet.Ni/ tBin *1e3;
    
    stairs(hAxe(3), tEdge(1:end-1), bSpkE, cSpecE(iter))
    stairs(hAxe(3), tEdge(1:end-1), bSpkI, cSpecI(iter))
end

% Plot the changed spike
if IdxSpk(IdxSpkDel) <= parsNet.Ne
    plot(hAxe(1), parsNet.T, parsNet.PrefStim(IdxSpk(IdxSpkDel)), 'co')
else
    plot(hAxe(2), parsNet.T, IdxSpk(IdxSpkDel) - parsNet.Ne, 'co')
end

axes(hAxe(1))
ylim(parsNet.PrefStim([1,end]))
set(gca, 'ytick', parsNet.PrefStim([1,end/2,end]), 'yticklabel', parsNet.Width*[-1,0,1])
set(gca, 'xtick', [])
ylabel('Neuron index')
axes(hAxe(2))
ylim([0, parsNet.Ni])
set(gca, 'xtick', [])
axes(hAxe(3))
xlim([0, parsNet.T])
xlabel('Time (ms)')
ylabel('Rate (Hz)')
ylim([0, max(bSpkI)*1.2])
linkaxes(hAxe, 'x')

% Specify the size of time bin
axes(hAxe(3))
axisRange = axis;
text(axisRange(2), axisRange(4), sprintf('tBin=%dms', tBin), ...
    'horizontalalignment', 'right', 'verticalalignment', 'top');
