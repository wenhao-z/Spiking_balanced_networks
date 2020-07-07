function rate = nlFunc_MeanField(mu, sigma, vth, tau, parsNet)
% The mean firing rate given the mean and flucatuation of synaptic input
% derived from mean field approximation
% Reference: Brunel, J. Comp. Neurosci., 2000, Equation 21.

% mu: the mean of synaptic input
% sigma: the std of synaptic input
% vth: the firing threshold 
% vre: reset potential after firing a spike
% tau: membrane time constant

% Rate: with unit 1/ms, i.e., kHz

% Author: Wen-Hao Zhang, Mar 12, 2019
% University of Pittsburgh
% wenhao.zhang@pitt.edu


% An approximation to overcome the zero initial condition
sigma(sigma==0) = sqrt(mu(sigma==0) * parsNet.FanoFactorInput); 

yth = (vth - mu) ./ sigma;
yre = (parsNet.vre - mu) ./ sigma;

Intgrand = @(u) (exp(u.*(-u + 2*yth)) - exp(u.*(-u + 2*yre)))./u;
% Note that when u==2*yth, an Inf value may appear

du = 0.01;
ugrid = du: du: 100;

rate = [2*(yth-yre), Intgrand(ugrid)]; 
% 2*(yth-yre) is the output of Intgrand when u=0, but the function
% handle is not able to deal with this singular condition

rate = sum(rate) * du;
rate = rate*tau + parsNet.refrac;
rate = 1./rate; % rau and refrac have unit of ms.

end