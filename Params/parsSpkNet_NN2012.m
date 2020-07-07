% Set the default network Parameters according to 
% Kumar & Doiron, Nat. Neurosci., 2012
% The parameters are packed into a struct

% Wen-Hao Zhang
% Math Department, University of Pittsburgh
% Jan-28, 2019

% ------------------------------------------------------------
parsNet.Ncells = 5000;
parsNet.Ne = 4000;
parsNet.Ni = 1000;
parsNet.Width = 180; % The width of feature space is 2*width, unit: degree
parsNet.TuneWidth = 40; % 40 degrees

PrefStim = linspace(-parsNet.Width, parsNet.Width, parsNet.Ne+1);
parsNet.PrefStim = PrefStim(2:end)';
clear PrefStim

% ------------------------------------------------------------
% Connection probabilities
parsNet.pee = 0.2;
parsNet.pie = 0.5;
parsNet.pei = 0.5;
parsNet.pii = 0.5;

% Scaled connection strength
parsNet.jee_scale = 10;
parsNet.jie_scale = 4;
parsNet.jei_scale = -16 *1.2;
parsNet.jii_scale = -16;
parsNet.ratio_jeestruct = 0.2; % the ratio of all structured E-E conns over all E-E conns
 
parsNet.gxe = 1; % gain of xe connections. x means any neuron type
parsNet.gxi = 1; % gain of xi connections. x means any neuron type

% parsNet.mue = 0*ones(1,2); % The range of resting potential for E neurons, [1.1, 1.2]
% parsNet.mui = 0*ones(1,2); % The range of resting potential for I neurons, [1, 1.05]

% ------------------------------------------------------------
% Feedforward inputs
parsNet.T = 2e3; %simulation time (ms)
parsNet.AmplIff = 0.1;
parsNet.Posi = 0; % The location of feedforward inputs in feature space. Unit: deg

parsNet.stimstart = 0; % (ms)
parsNet.stimend = 0; % (ms)

% Background inputs
parsNet.IBkg = [0.07; 0.1]; % The background inputs applied to E and I neurons respectively
parsNet.FanoFactorInput = 1;

% ------------------------------------------------------------
parsNet.maxrate = 100; %(Hz) maximum average firing rate.
% if the average firing rate across the simulation for any neuron exceeds
% this value, some of that neuron's spikes will not be saved

% Set the random number seed
parsNet.rngWeight = rng('shuffle');
parsNet.rngInput = rng('shuffle');
% parsNet.seedSim = sum(clock)*100;
% parsNet.rsWeight = RandStream.create('mt19937ar', 'Seed', parsNet.seedWeight);


%% Intrinsic parametrs of single neuorns
% Do not CHANGE these parameters !!
parsNet.taue = 15; %membrane time constant for exc. neurons (ms)
parsNet.taui = 10;

% parsNet.EL      = 0; This reminds us that the resting potential is 0mV, because
%              it doesn appear in the iterative dynamics
parsNet.vre     = 0; %reset voltage
parsNet.threshe = 1; %threshold for exc. neurons
parsNet.threshi = 1;

parsNet.dt = 0.1; %simulation timestep (ms)
parsNet.refrac = 5; %refractory period (ms)

%synaptic time constants (ms)
parsNet.tauerise    = 1;
parsNet.tauedecay   = 3;
parsNet.tauirise    = 1;
parsNet.tauidecay   = 2;

%%

parsNet = orderfields(parsNet);