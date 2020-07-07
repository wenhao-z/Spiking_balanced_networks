function weights = genNetConnMat(parsNet)
% Generate the connection matrix for spiking network simulation
% Author: Wen-Hao Zhang
% wenhao.zhang@pitt.edu
% University of Pittsburgh
% Feb 5, 2019

% The connection weight has two parts
% 1) sparse random 
% 2) structured connection on a ring between E-E neurons
% It is probabilistic to determine whether two neurons are connected or
% not. Once two neurons are connected, their connection weight will be set
% according to the weight jxx. E-I, I-I, I-E connections are not
% structured, but E-E connections have a ring structure whose connection
% strength depends on the distance between two neurons on the ring.

Ne  = parsNet.Ne;
Ni  = parsNet.Ni;
Jmat = reshape(parsNet.Jmat, 2, 2);
connProb = reshape(parsNet.connProb, 2, 2);

% Load structured connection parameter
TuneWidth  = parsNet.TuneWidth;
Width      = parsNet.Width;

% E-E connection matrix
W = parsNet.PrefStim - parsNet.PrefStim(1);
W = angle(exp(W*pi*1i/Width))* Width/pi;
W = exp(-W.^2/(2*TuneWidth^2))/(sqrt(2*pi)*TuneWidth);
W = gallery('circul', W);

ratio_jeestruct = parsNet.ratio_jeestruct;
W = Jmat(1,1) *(ratio_jeestruct*2*Width* W + 1 - ratio_jeestruct);

% set the state of random number generator
rng(parsNet.rngWeight); 

% Full connection matrix
weights = [W.* (rand(Ne,Ne) < connProb(1,1)), Jmat(1,2)*(rand(Ne,Ni) < connProb(1,2)); ...
    Jmat(2,1)*(rand(Ni,Ne) < connProb(2,1)), Jmat(2,2)*(rand(Ni,Ni) < connProb(2,2))];

% No self connections
weights(1:1+size(weights,1):end) = 0;

end
