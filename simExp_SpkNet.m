% Simulate the experiments to generate the results of decentralized network
% model consisting of congruent and opposite neurons.

% Author: Wen-Hao Zhang, Feb. 7, 2019
% wenhao.zhang@pitt.edu
% @ University of Pittsburgh

setWorkPath;
addpath(fullfile(Path_RootDir, 'simExp'));
addpath(fullfile(Path_RootDir, 'SpikeAnalysisKit'));

%% Load DEFAULT parameters
parsSpkNet;

%% Perform virtual experiments
flagExp = 6;
% 1. Demo the membrane potential trajectory and rastergram
% 2. Demo the chaotic state of network through comparing the network
%    activity with and without deleting a spike.
% 3. Make a input-firing rate curve of the network 
% 4. Test the firing rate and its distribution with connection strength in
%    a homogeneous network
% 5. Linear stability analysis of network
% 6. Demo the spiking activity of network with ring structure
% 7. Analyze the fluctuations of network activity on ring manifold and test
%    its dependence with connection strength

switch flagExp
    case 1
        demoNet_PotentialSpks;
    case 2
        demoChaoticState;
    case 3
        getNetIOCurve;
    case 4
        demoNetSpks_ConnStr;
    case 5
        ansLinearStability;
    case 6
        demoNet_RingConn;
    case 7
        ansNet_RingManifold;
end
