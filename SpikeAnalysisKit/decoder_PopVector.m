function [DecodeRes, DecodePars] = decoder_PopVector(tSpk, parsNet, varargin)
% % Using population vector to decode the network activity on a ring
% manifold, and calculate the statistics of the decoded results.

% Wen-Hao Zhang
% Math Department, University of Pittsburgh
% wenhao.zhang@pitt.edu
% Feb 8, 2019

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

%% Population vector decoding
neuronEdge = [1: parsNet.Ne, parsNet.Ne + 0.5];
if ~exist('tStat', 'var')
    tStat = 100; % unit: ms
end
if ~exist('tBin', 'var')
    tBin = 5; % unit: ms
end
tEdge = tStat: tBin : parsNet.T;

bSpk = histcounts2(tSpk(1,:)', tSpk(2,:)', neuronEdge, tEdge);

cirPos = exp(1i * parsNet.PrefStim/ parsNet.Width * pi);

BumpPos = mean(bsxfun(@times, cirPos, bSpk), 1); % average over all neurons
% Note that angle(0)=0, but here 0 comes from no spikes at all.
BumpPos = BumpPos ./ abs(BumpPos); % normalize
if exist('bOutBumpPos', 'var') && bOutBumpPos
    DecodeRes.BumpPos = angle(BumpPos)* parsNet.Width/pi; % angular value
end

BumpPos(isnan(BumpPos)) = []; % Remove the time bin without spikes

% Statistics of bump position
meanBumpPos = mean(BumpPos, 2);
mrlBumpPos = abs(meanBumpPos); % mean resultant length, average over time
concBumpPos = mrl2Kappa(mrlBumpPos); % concentration parameter

meanBumpPos = angle(meanBumpPos)* parsNet.Width/pi; % angular value
BumpPos = angle(BumpPos)* parsNet.Width/pi; % angular value

% Variance
devBumpPos = bsxfun(@minus, BumpPos, meanBumpPos);
devBumpPos(devBumpPos>parsNet.Width) = devBumpPos(devBumpPos>parsNet.Width) - 2*parsNet.Width;
devBumpPos(devBumpPos<parsNet.Width) = devBumpPos(devBumpPos<parsNet.Width) + 2*parsNet.Width;
varBumpPos = cov(devBumpPos');


%% Fold output arguments into a struct
% Output the parameters for statistics
DecodePars.tBin     = tBin;
DecodePars.tStat    = tStat;

DecodeRes.meanBumpPos   = meanBumpPos;
DecodeRes.concBumpPos   = concBumpPos;
DecodeRes.mrlBumpPos    = mrlBumpPos;
DecodeRes.varBumpPos    = varBumpPos;

end