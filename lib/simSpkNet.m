function outSet = simSpkNet(parsNet, outArgsOpt, varargin)
% This code simulates a spiking neural network with E and I neurons
% Author: Wen-Hao Zhang
% Math Department, University of Pittsburgh
% wenhao.zhang@pitt.edu
% Jan 30, 2019

% fprintf("setting up parameters\n")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the possible parameters from varargin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PS: I typically use this to convey v and bSpk, which fully characterizes
% the network state. 
if mod(size(varargin,2) , 2) == 1 % odd number input
    error('The varargin input number is wrong!')
else
    for iter = 1: round(size(varargin,2)/2)
        eval([varargin{2*iter-1} '= varargin{2*iter};']);
    end
end
clear iter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Release some common parameters from parsNet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ncells 	= parsNet.Ncells;
Ne 		= parsNet.Ne;
Ni 		= parsNet.Ni;
T		= parsNet.T;

vre     = parsNet.vre; %reset voltage
dt      = parsNet.dt; %simulation timestep (ms)
refrac  = parsNet.refrac; %refractory period (ms)

% Some auxiliary variables to speed up the simulation of PSP
alphaerise  = 1 - dt/parsNet.tauerise;
alphaedecay = 1 - dt/parsNet.tauedecay;
alphairise  = 1 - dt/parsNet.tauirise;
alphaidecay = 1 - dt/parsNet.tauidecay;

df_taue = parsNet.tauedecay - parsNet.tauerise;
df_taui = parsNet.tauidecay - parsNet.tauirise;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Connection matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
weights = genNetConnMat(parsNet);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the stimulus/input parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stimstart = parsNet.stimstart;
% stimend   = parsNet.stimend;

IBkg = parsNet.IBkg;
IBkg = [IBkg(1)*ones(Ne,1); IBkg(2)*ones(Ni,1)];

stdIBkg = sqrt(IBkg * parsNet.FanoFactorInput);
stdIBkg = [stdIBkg(1)*ones(Ne,1); stdIBkg(2)*ones(Ni,1)];

Iff = makeIff(parsNet);
Iff = [Iff; zeros(Ni,1)]; % divide by time constant tau to speed up iteration
% Iff = [Iff/parsNet.taue; zeros(Ni,1)]; % divide by time constant tau to speed up iteration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize arrays for network simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thresh = [parsNet.threshe*ones(Ne,1); parsNet.threshi*ones(Ni,1)];
tau    = [parsNet.taue*ones(Ne,1);    parsNet.taui*ones(Ni,1)];
alphav = 1-dt./tau; % Auxiliary variable for updating membrane potential v

maxTimes = round(parsNet.maxrate*T/1000);
tSpk     = zeros(2, Ncells*maxTimes);

xerise  = zeros(Ncells,1); % auxiliary variables for E/I currents (difference of exponentials)
xedecay = zeros(Ncells,1);
xirise  = zeros(Ncells,1);
xidecay = zeros(Ncells,1);

if ~exist('v', 'var')
    v = rand(Ncells,1); %membrane voltage
end
if ~exist('bSpk', 'var')
    bSpk = false(Ncells,1);
end

lastSpike = -100*ones(Ncells,1); %time of last spike
Nsteps = round(T/dt);

nSpkAll = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf('Starting simulation')
rng(parsNet.rngInput); % set the state of random number generator
for ti = 1:Nsteps
    t = dt*ti;
    
    % Recurrent input (each PSP satisfies an alpha function)
    IrecE = sum(weights(:, bSpk(1:Ne)), 2);
    IrecI = sum(weights(:, [false(Ne,1); bSpk((Ne+1):end)]), 2);
    
    xerise 	= alphaerise  * xerise  + IrecE;
    xedecay = alphaedecay * xedecay + IrecE;
    xirise 	= alphairise  * xirise  + IrecI;
    xidecay = alphaidecay * xidecay + IrecI;
    
    Isyn = (xedecay - xerise)/df_taue + (xidecay - xirise)/df_taui;
        
    % Apply feedforward inputs and background inputs
    Isyn = Isyn + Iff + IBkg + stdIBkg.*randn(Ncells,1)/sqrt(dt);
    
    %     if (t > stimstart) && (t < stimend)
    %         Isyn(1:Nstim) = Isyn(1:Nstim) + stimstr;
    %     end
    
    % find the neurons not in refractory period
    idxCell = (t> (lastSpike+refrac));
    v(idxCell) = alphav(idxCell).*v(idxCell) + Isyn(idxCell)*dt;
    
    bSpk = (v > thresh); % binary spike pattern
    v(bSpk) = vre;
    lastSpike(bSpk) = t;
    
    % Record the spike timing
    nSpkNow = length(find(bSpk));
    tSpk(1, nSpkAll+ (1:nSpkNow)) = find(bSpk); % Index of spiking neurons
    tSpk(2, nSpkAll+ (1:nSpkNow)) = t; % Spike timing
    nSpkAll = nSpkAll + nSpkNow;
    
end %end loop over time
% fprintf('\r')

%% Fold the outputs into a struct
% Delete zeros in tSpk
tSpk(:, nSpkAll+1:end) = [];
outSet.tSpk = tSpk;

% outSet.nSpk = accumarray(tSpk(1,:)', 1);
% outSet.weights = sparse(weights);

if exist('outArgsOpt', 'var')
    if isfield(outArgsOpt, 'v')
        outSet.v = v;
    end
    if isfield(outArgsOpt, 'bSpk')
        outSet.bSpk = bSpk;
    end
end

end


function Iff = makeIff(parsNet)
% Generate the mean feedforward input whose profile is a Gaussian bump 
Width     = parsNet.Width;
TuneWidth = parsNet.TuneWidth;
POSI      = parsNet.Posi;

POSI = POSI- parsNet.PrefStim;

% Add periodic condition
POSI = angle(exp(1i*POSI * pi/Width)) * Width/pi;
Iff = exp(-(POSI).^2/ (4*TuneWidth^2)); % [N, 1]
Iff = Iff * parsNet.AmplIff /(2*sqrt(pi)*TuneWidth); % Normalized. The sum of Iff over neurons is AmplIff
Iff = Iff * 2 * Width; % Normalized. The sum of Iff is Ne * AmplIff

end