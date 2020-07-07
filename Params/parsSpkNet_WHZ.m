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
parsNet.pee = 0.2;
parsNet.pei = 0.2;
parsNet.pie = 0.2;
parsNet.pii = 0.2;

% Scaled connection strength
parsNet.jxe = 2; % E synaptic weight
parsNet.ratiojie = 5; % The ratio between I synapse over E synapse, jxi/jxe. (absoluate value)
parsNet.ratio_jeestruct = 0; % the ratio of all structured E-E conns over all E-E conns

sqrtN = sqrt(parsNet.Ncells);
parsNet.jee = parsNet.jxe/sqrtN;
parsNet.jie = parsNet.jxe/sqrtN;
parsNet.jei = -parsNet.ratiojie * parsNet.jxe/sqrtN * 1.1;
parsNet.jii = -parsNet.ratiojie * parsNet.jxe/sqrtN;

% ------------------------------------------------------------
% Feedforward inputs
parsNet.T = 2e3; %simulation time (ms)
parsNet.AmplIff = 7.5e-3; % Peak intensity of bump feedforward input.
%                           Unit: mV. => IBkg will be divided by tau in iterative dynamics
parsNet.Posi = 0; % The location of feedforward inputs in feature space. Unit: deg

parsNet.stimstart = 0; % (ms)
parsNet.stimend = 0; % (ms)

% Background inputs
parsNet.IBkg_scale = [1;1]; % The background inputs applied to E and I neurons respectively
%                             IBkg is scaled by thresh/tau_m with unit kHz
parsNet.IBkg = parsNet.IBkg_scale .* ...
    [parsNet.threshe/parsNet.taue; parsNet.threshi/parsNet.taui];

parsNet.FanoFactorInput = 0.01; % This corresponds to the feedforward connection weight


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

