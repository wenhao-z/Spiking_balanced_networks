function ISIRes = intSpkIntval(tSpk, parsNet)
% Calculate the inter-spike interval (ISI) of spikes emitted by the same
% type of neuron in a HOMOGENEOUS network model

% Wen-Hao Zhang
% Math Department, University of Pittsburgh
% wenhao.zhang@pitt.edu
% Feb 5, 2019

%           Input structure
% tSpk: 1st row is the index of neurons
%       2nd row is the spike timing (unit: ms)

% Change the structure of tSpk from
%  a [2, nSpkAll] array        into
%  a cell[Ncells, 1] with each element record the spike times emitted by a neuron
if isempty(tSpk)
    ISIRes = struct('histISIE', [], 'hist_ISIO', [], ...
        'edgeISIE', [], 'edgeISII', [], ...
        'CV_ISIE', [], 'CV_ISII', []);
    
    return
end
tSpkCell = accumarray(tSpk(1,:)', tSpk(2,:)', [parsNet.Ncells, 1], @(x) {sort(x)});
tISI = cellfun(@diff, tSpkCell, 'uniformoutput', 0);
clear tSpkCell

tISIE = cell2mat(tISI(1:parsNet.Ne))';
tISII = cell2mat(tISI(parsNet.Ne+1:end))';

%% Histogram of Inter-spike intervals of E and I neurons
nBins = 50;
[histISIE, edgeISIE] = histcounts(tISIE, nBins, 'Normalization', 'probability');
[histISII, edgeISII] = histcounts(tISII, nBins, 'Normalization', 'probability');

%% CV (coefficient of variation) of Inter-Spike interval
% CV = standard deviation / mean
CV_ISIE = std(tISIE)/mean(tISIE);
CV_ISII = std(tISII)/mean(tISII);

%% Fold results into an output struct
ISIRes.histISIE = histISIE;
ISIRes.hist_ISIO = histISII;
ISIRes.edgeISIE = edgeISIE;
ISIRes.edgeISII = edgeISII;

ISIRes.CV_ISIE = CV_ISIE;
ISIRes.CV_ISII = CV_ISII;