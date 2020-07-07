% The parameters of spiking neural network
% Wen-Hao Zhang, Oct-24, 2016
% wenhaoz1@andrew.cmu.edu
% @Carnegie Mellon University

% Ref: Ashok Litwin-Kunar, Nat. Comms. 2014, Formation and maintenance of
% neuronal assemblies through synaptic plasticity

%% Membrane dynamics
NetPars.tauE        = 20; % E membrane time constant, ms   taue
NetPars.tauI        = 20; % I membrane time constant, ms   taui
NetPars.restV_E     = -70; % E resting potential, mV  vleake
NetPars.restV_I     = -62; % I resting potential, mV  vleaki
NetPars.deltaTheta  = 2; % Exponential IF slope parameter, mV deltathe
NetPars.C           = 300; % Capacitance, pF

NetPars.revV_E      = 0; % E synapse reversal potential, mV  erev
NetPars.revV_I      = -75; % I synapse reversal potntial, mV  irev
NetPars.thetaV      = -52; % initial spike voltage threshold, mV vth0
NetPars.ATheta      = 10; % Increase in threshold post spike, mV ath
NetPars.tauTheta	= 30; % Threshold decay timescale, ms tauth
NetPars.resetV      = -60; % Reset potential after firing a spike, mV vre
NetPars.tauRefrac   = 1; % Absolute refractory period, ms taurefrac
NetPars.aIadapt     = 4; % Adaptation parameter a, nS  aw_adapt
NetPars.bIadapt     = .805; % Adaptation parameter b, pA
NetPars.tauIadapt   = 150; % Adaptation timescale  tauw_adapt

%% Connectivity
NetPars.Ne          = 4000; % Number of E neurons
NetPars.Ni          = 1000; % Number of I neurons
NetPars.Ncells      = NetPars.Ne + NetPars.Ni;
NetPars.Jee0        = 2.86; % Initial E-to-E strength
NetPars.Jei0        = 48.7; % Initial I-to-E strength
NetPars.Jie         = 1.27; % E-to-I strength (not plastic)
NetPars.Jii         = 16.2; % I-to-I strength (not plastic)
NetPars.connProb    = 0.2;
    
NetPars.tauRiseE    = 1; % E synapse rise time tauerise, ms
NetPars.tauDecayE   = 6; % E synapse decay time tauedecay, ms
NetPars.tauRiseI    = .5; % I synapse rise time tauirise, ms
NetPars.tauDecayI   = 2; % I synapse decay time tauidecay, ms

NetPars.rateEF      = 4.5; % Feedforward input rate to e (khz)  rex
NetPars.rateIF      = 2.25; % Feedforward input rate to i (khz)  rix

NetPars.JeeMin      = 1.78; % Minimum ee strength
NetPars.JeeMax      = 21.4; % Maximum ee strength

NetPars.JeiMin      = 48.7; % Minimum ei strength
NetPars.JeiMax      = 243; % Maximum ei strength

NetPars.Jef         = 1.78; % Feedforward to E strength Jex  
NetPars.Jif         = 1.27; % Feedforward to i strength Jix

%% Voltage based STDP (modifying E-to-E connections)
NetPars.aLTD        = .0008; % LTD strength  altd
NetPars.aLTP        = .0014; % LTP strength  altp
NetPars.thetaLTD    = -70; % LTD voltage threshold  thetaltd
NetPars.thetaLTP    = -49; % LTP voltage threshold  thetaltp
NetPars.tauu        = 10; % Timescale for u variable
NetPars.tauv        = 7; % Timescale for v variable
NetPars.taux        = 15; % Timescale for x variable

%% Inhibitory STDP (modifying I-to-E connections)
NetPars.tauy        = 20; % Width of istdp curve
NetPars.eta         = 1; % iSTDP learning rate
NetPars.r0          = .003; % Target rate (khz)

%% Populations
NetPars.Npop        = 20; %number of assemblies
NetPars.pmembership = .05; %probability of belonging to any assembly
NetPars.Nmaxmembers = 300; %maximum number of neurons in a population (to set size of matrix)

%% Simulation
NetPars.dt          = 0.1; % Integration timestep, ms
NetPars.tLen        = 2e3; % Simulation time, ms
NetPars.Nskip       = 1000; % How often (in number of timesteps) to save w_in
NetPars.Vpeak       = 20; % Cutoff for voltage.  when crossed, record a spike and reset, mV
NetPars.dtnormalize = 20; % How often to normalize rows of ee weights, ms?
NetPars.stdpdelay   = 1000; % Time before stdp is activated, to allow transients to die out
NetPars.Nspikes     = 100; % Maximum number of spikes to record per neuron