function [rate, Isyn, VarIsyn] = getFiringRateTheory(parsNet)
% Predict the mean firing rate of a spiking neural network

% Author: Wen-Hao Zhang
% wenhao.zhang@pitt.edu
% University of Pittsburgh
% Feb 9, 2019

% Reference: Brunel, J. Comp. Neurosci., 2000

% Effective connection matrix in calculating MEAN input after mean-field approximation
% Jxy,eff = sqrt(N) * Pxy * Cy * Jxy,
% where N is the amount of neurons in the network
%       Pxy is the connection probability from neural population b to a
%       Cy is the fraction of type-y neuron, i.e., Cy = Ny/N
%       Jxy is the mean synaptic weight from neural population b to a

% Effective connection matrix in calculating VARIANCE of input after mean-field approximation
% Sigxy,eff = Pxy * Cy * Jxy^2,


% Rate: unit Hz
% Isyn: unit is 1/ms, i.e., kHz
% VarIsyn: unit is 1/ms^2

Cmat = [parsNet.Ne, parsNet.Ni];
Jmat = reshape(parsNet.Jmat, 2, 2);
CP = reshape(parsNet.connProb,2,2) .* Cmat(:)'; % number of connections. Ensure that Cmat is a row vector

JeffMat = Jmat .* CP;
VarMat = Jmat.^2 .* CP;

% Feedforward input
Iff = parsNet.IBkg; % Note the unit should be 1/ms. This is consistent with tau with ms

% Initialize the firing rate, one for E neurons and another for I neurons
rate = zeros(2,1);
rateNew = zeros(2,1);

% Iteratively solve for rate
drate = 10;
dt = 0.05;
eps = 1e-9;
maxIter = 1e5;

nIter = 1;
while sum(abs(drate)) > eps && nIter < maxIter
    % Note that during iteration, rate has unit of ms.
    Isyn    = JeffMat * rate + Iff;
    VarIsyn = VarMat * rate + Iff * parsNet.FanoFactorInput;
    
    Isyn = Isyn .* [parsNet.taue; parsNet.taui];
    VarIsyn = VarIsyn .* [parsNet.taue; parsNet.taui];
    
    % Input-firing rate nonlinearity after mean-field approximation
    rateNew(1) = nlFunc_MeanField(Isyn(1), sqrt(VarIsyn(1)), parsNet.threshe, parsNet.taue, parsNet);
    rateNew(2) = nlFunc_MeanField(Isyn(2), sqrt(VarIsyn(2)), parsNet.threshi, parsNet.taui, parsNet);
    
    drate = rateNew - rate;
    rate = rate + drate*dt;
    %     drate = rateNew - rate(:,end);
    %     rate = [rate, rate(:,end) + drate*dt];
    
    nIter = nIter + 1;
end
rate = rate * 1e3; % convert the unit from 1/ms to 1/s

end