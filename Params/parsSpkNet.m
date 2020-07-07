% Set the default network Parameters
% The parameters are packed into a Dictionary

% Wen-Hao Zhang
% Math Department, University of Pittsburgh
% Jan-28, 2019

% ------------------------------------------------------------
parsNet.Ncells  = 5000;
parsNet.Ne      = 4000;
parsNet.Ni      = 1000;
parsNet.Width   = 180; % The width of feature space is 2*width, unit: degree
parsNet.TuneWidth = 40; % 40 degrees

PrefStim = linspace(-parsNet.Width, parsNet.Width, parsNet.Ne+1);
parsNet.PrefStim = PrefStim(2:end)';
clear PrefStim

% ------------------------------------------------------------
% Connection probabilities
% connProb = [pee, pei; pie, pii]
parsNet.connProb = [0.2, 0.2; 0.2, 0.2];
parsNet.connProb = parsNet.connProb(:); % Reshape it into a column vector in consistent with paramGrid.m

% Scaled connection strength
parsNet.jxe = 15; % E synaptic weight 
% parsNet.jxi = 16; % I synaptic weight (absoluate value)
parsNet.ratiojie = 5; % The ratio between I synapse over E synapse, jxi/jxe. (absoluate value)
parsNet.ratio_jeestruct = 0; % the ratio of all structured E-E conns over all E-E conns

% Uncomment following 4 lines could give more freedom of connection stregnth.
% Following is the parameter from Ashok, Nat. Neurosci., 2012
% parsNet.jee_scale = 10;
% parsNet.jie_scale = 4;
% parsNet.jei_scale = -16*1.2;
% parsNet.jii_scale = -16;

% parsNet.mue = 0*ones(1,2); % The range of resting potential for E neurons, [1.1, 1.2]
% parsNet.mui = 0*ones(1,2); % The range of resting potential for I neurons, [1, 1.05]

% ------------------------------------------------------------
% Feedforward inputs
parsNet.T = 2e3; %simulation time (ms)
parsNet.AmplIff = 7.5e-3; % Peak intensity of bump feedforward input.
%                           Unit: mV. => IBkg will be divided by tau in iterative dynamics
parsNet.Posi = 0; % The location of feedforward inputs in feature space. Unit: deg

parsNet.stimstart = 0; % (ms)
parsNet.stimend = 0; % (ms)

% Background inputs
parsNet.IBkg_scale = [1.1; 1.05]; % The background inputs applied to E and I neurons respectively
%                             IBkg is scaled by thresh/tau_m with unit kHz

parsNet.FanoFactorInput = 0;

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
%                      it doesn't appear in the iterative dynamics
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