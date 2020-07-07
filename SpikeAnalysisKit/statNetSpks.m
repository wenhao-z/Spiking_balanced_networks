function [NetStatRes, NetStatPars] = statNetSpks(tSpk, parsNet, listAnalysis, varargin)
% % Using population vector to decode the network activity on a ring
% manifold, and calculate the statistics of the decoded results.

% Wen-Hao Zhang
% Math Department, University of Pittsburgh
% wenhao.zhang@pitt.edu
% Feb 8, 2019

%% --------------------------------------------------------
% Get possible parameters from varargin
% Use this to convey the parameter tBin, and bInit
if mod(size(varargin,2) , 2) == 1 % odd number input
    error('The varargin input number is wrong!')
else
    for iter = 1: round(size(varargin,2)/2)
        eval([varargin{2*iter-1} '= varargin{2*iter};']);
    end
end
clear iter

if ~exist('bInit', 'var')
    bInit = 0;
end

%% --------------------------------------------------------
% Initialize the struct of NetStatRes
if bInit == 1
    for var = listAnalysis
        switch var{1}
            case 'PopVector'
                % Use population vector to decode the network's bump position
                % on ring manifold, and calculate the statistics of bump
                % position
                NetStatRes.meanBumpPos  = [];
                NetStatRes.concBumpPos  = [];
                NetStatRes.mrlBumpPos   = [];
                NetStatRes.varBumpPos   = [];
            case 'IntSpkIntval'
                % Calculate the inter-spike interval of neurons with the same type
                NetStatRes.ISIRes = [];
            case 'Rate_AvgXTime'
                % Firing rate of all neurons averaged over TIME
                NetStatRes.fireRate_neuron = [];
            case 'Rate_AvgXNeuron'
                % Firing rate of all neurons averaged over NEURON
                NetStatRes.fireRate_time = [];
                NetStatRes.tEdge_Rate    = [];
            case 'fanoFactor_nSpk'
                % Fano factor of spike count of all neurons
                NetStatRes.fanoFactor_nSpk = [];
        end
    end
    
    return
end

%% --------------------------------------------------------
% Perform the analysis of network spikes
if ~exist('tStat', 'var')
    tStat = 100; % unit: ms
end
if ~exist('tBin', 'var')
    tBin = 5; % unit: ms
end

for var = listAnalysis
    switch var{1}
        case 'PopVector'
            % Use population vector to decode the network's bump position
            % on ring manifold, and calculate the statistics of bump
            % position
            [DecodeRes, DecodePars] = decoder_PopVector(tSpk, parsNet, ...
                'tStat', tStat, 'tBin', tBin);
            for varName = fieldnames(DecodeRes)'
                NetStatRes.(varName{1}) = DecodeRes.(varName{1});
            end
            NetStatPars.DecodePars = DecodePars;
            
        case 'IntSpkIntval'
            % Calculate the inter-spike interval of neurons with the same type
            ISIRes = intSpkIntval(tSpk, parsNet);
            NetStatRes.ISIRes = ISIRes;
            
        case 'Rate_AvgXTime'
            % Firing rate of all neurons averaged over TIME
            NetStatRes.fireRate_neuron = meanFiringRate(tSpk, parsNet, 'time', ...
                'tStat', tStat);
            
        case 'Rate_AvgXNeuron'
            % Firing rate of all neurons averaged over NEURON
            [fireRate, tEdge_Rate] = meanFiringRate(tSpk, parsNet, 'neuron', 'tBin', 5);
            NetStatRes.fireRate_time = fireRate;
            NetStatRes.tEdge_Rate = tEdge_Rate;
        case 'fanoFactor_nSpk'
            % Fano factor of spike count of all neurons
            NetStatRes.fanoFactor_nSpk = fanoFactorNSpk(tSpk, parsNet);
            
    end
end

end