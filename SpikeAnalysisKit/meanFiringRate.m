function [fireRate, tEdge] = meanFiringRate(tSpk, parsNet, dim_avg, varargin)
% Calculate the mean firing rate of every neurons in a network model

% Wen-Hao Zhang
% Math Department, University of Pittsburgh
% wenhao.zhang@pitt.edu
% Feb 8, 2019

%           Input structure
% tSpk: 1st row is the index of neurons
%       2nd row is the spike timing (unit: ms)
% dim_avg: the dimension over which the average is taking along. 
%          There are two choices:
%          1): 'time' (default): average over time. 
%              The result of the mean firing rate of all neurons
%          2): 'neuron': averge over neurons of the same type. 
%              The result is the mean firing rate over time.

% --------------------------------------------------------
% Get possible parameters from varargin
% Use this to convey the parameter tBin
if mod(size(varargin,2) , 2) == 1 % odd number input
    error('The varargin input number is wrong!')
else
    for iter = 1: round(size(varargin,2)/2)
        eval([varargin{2*iter-1} '= varargin{2*iter};']);
    end
end
clear iter
% --------------------------------------------------------

if ~exist('dim_avg', 'var')
   dim_avg = 'time';
end


switch dim_avg
    case 'time'
        % Mean firing rate of all neurons
        if ~exist('tStat', 'var')
            tStat = 0; % unit: ms
        end
        
        % Remove the spikes before tStat, which is transient response
        tSpk(:, tSpk(2,:)<tStat) = [];
        tSpk(2,:) = tSpk(2,:) - tStat;
        
        nSpk = accumarray(tSpk(1,:)', 1);
        if length(nSpk) < parsNet.Ncells
            % When the neurons with large index number don't fire any spike, the
            % lenght of nSpk could be smaller than Ncells
            nSpk = [nSpk; zeros(parsNet.Ncells-length(nSpk),1)];
        end
        fireRate = nSpk/(parsNet.T - tStat)*1e3; % Unit: Hz
        
    case 'neuron'
        if ~exist('tBin', 'var')
            tBin = 5; % unit: ms
        end
        tEdge = [0: tBin : parsNet.T, parsNet.T + tBin/2];
        
        neuronEdge = [1, parsNet.Ne + 0.5];
        rateE = histcounts2(tSpk(1,:)', tSpk(2,:)', neuronEdge, tEdge);
        rateE = rateE/ parsNet.Ne/ tBin *1e3;
        
        neuronEdge = [parsNet.Ne+1, parsNet.Ncells+0.5];
        rateI = histcounts2(tSpk(1,:)', tSpk(2,:)', neuronEdge, tEdge);
        rateI = rateI/ parsNet.Ni/ tBin *1e3;
        
        fireRate = [rateE; rateI];        
end

