function fanoFactor_nSpk = fanoFactorNSpk(tSpk, parsNet, varargin)
% Calculate the Fano factor of spike count emitted during different segment.

% Wen-Hao Zhang
% Math Department, University of Pittsburgh
% wenhao.zhang@pitt.edu
% Feb 8, 2019

%           Input structure
% tSpk: 1st row is the index of neurons
%       2nd row is the spike timing (unit: ms)

% --------------------------------------------------------
% Get possible parameters from varargin
% Use this to convey the parameter tSeg
if mod(size(varargin,2) , 2) == 1 % odd number input
    error('The varargin input number is wrong!')
else
    for iter = 1: round(size(varargin,2)/2)
        eval([varargin{2*iter-1} '= varargin{2*iter};']);
    end
end
clear iter
% --------------------------------------------------------


neuronEdge = [1: parsNet.Ncells, parsNet.Ncells+0.5];
tSeg = 0.5e3; % duration of each segment, unit: ms
tEdge = tSeg: tSeg: parsNet.T; % Omit the 1st segment from analysis to rule out the initial condition effect

nSpk = histcounts2(tSpk(1,:)', tSpk(2,:)', neuronEdge, tEdge);

rateAvg = mean(nSpk,2); 
fanoFactor_nSpk = var(nSpk, 0, 2) ./ rateAvg;
fanoFactor_nSpk(rateAvg==0) = 0;
