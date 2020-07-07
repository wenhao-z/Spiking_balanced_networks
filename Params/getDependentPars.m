function parsNet = getDependentPars(parsNet)
% Calculate dependent network parameters
% Wen-Hao Zhang,
% wenhao.zhang@pitt.edu
% University of Pittsburgh
% Feb 5, 2019

if isfield(parsNet, 'jxi') && isfield(parsNet, 'ratiojie')
    error('Specify either jxi or ratiojie, but not both!')
end

% Jmat_scale is the connection strength which is invariant with
% __time constant of membrane potential__ and __network size__
% Jmat = [Jee, Jei; Jie, Jii]
if isfield(parsNet, 'jxe') && isfield(parsNet, 'ratiojie')
    Jmat_scale = parsNet.jxe * repmat([1, -parsNet.ratiojie], 2, 1);
end
if isfield(parsNet, 'jxe') && isfield(parsNet, 'jxi')
    Jmat_scale = repmat([jxe, -jxi], 2, 1);
end
% Jmat_scale(1,1) = Jmat_scale(1,1) * 1.2;  % Jee
Jmat_scale(1,2) = Jmat_scale(1,2) * 1.2;  % Jei


% The actual synaptic weight used to generate the raw weight matrix
% Normalize the weight by the square root of neuron number, which is
% necessary condition to achieve chaotic balanced state
Jmat = Jmat_scale / sqrt(parsNet.Ncells)...
    ./ [parsNet.taue; parsNet.taui];

% The parameter equivalent to Ashok's Nat. Neurosci., 2010 paper.
% parsNet.Jmat = [0.0236, -0.1131; 0.0355 -0.1414];


% Reshape Jmat into column vectors, in order to be consistent with
% paramGrid.m, because it detects the vector with different columns and
% deem it as multiple values of parameters.
parsNet.Jmat_scale = Jmat_scale(:);
parsNet.Jmat = Jmat(:);

% ------------------------------------------------------------
% Feedforward input strength
% Input strength is scaled by the value
%    firing threshold / membrane time constant
% IBkg has unit of kHz, because tau has unit of ms
parsNet.IBkg = parsNet.IBkg_scale .* ...
    [parsNet.threshe/parsNet.taue; parsNet.threshi/parsNet.taui];

end