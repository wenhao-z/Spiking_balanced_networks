function [tSpk, W, VArray, bSpkFFAarray, IA] = simSpkNet(NetPars)
% Simulate the spiking neural network with Excitatory and inhibitory neurons
% Wen-Hao Zhang, Oct-24, 2016
% wenhaoz1@andrew.cmu.edu
% @Carnegie Mellon University

% OUTPUTS
% tSpk: timings of spikes, ms
% W:    connection matrix
% V:    membrane potential

% Unfold parameters from NetPars
for varName = fieldnames(NetPars)'
    eval([varName{1} '= NetPars.(varName{' num2str(1) '});']);
end

%% Initialization
nStep = round(tLen/dt);
V = resetV + (thetaV-resetV)*rand(NetPars.Ncells, 1);
% V = resetV * ones(NetPars.Ncells, 1);
if nargout > 2
    VArray = zeros(NetPars.Ncells, nStep); % Membrane potential
    VArray(:,1) = V;
end
if nargout >3
    bSpkFFAarray = zeros(NetPars.Ncells, nStep);
end

if nargout > 4
    IA = zeros(NetPars.Ne, nStep);
end

vth = thetaV * ones(NetPars.Ne, 1); % dynamic firing threshold for E neuron

tLastSpk = -10*ones(NetPars.Ncells, 1);
tSpk = zeros(2, 200*Ncells); % Spike timing
nSpkAll = 1; % number of all spikes of all neurons
bSpkNet = false(NetPars.Ncells, 1); % Binary response of last time step

% Feedforward inputs
rateExt = [rateEF * ones(NetPars.Ne, 1); ...
    rateIF * ones(NetPars.Ni, 1)] * dt; % unit: hz/bin
Jff = [Jef * ones(NetPars.Ne, 1); Jif * ones(NetPars.Ni, 1)];

% Adaptation current, only for E neurons
Iadapt = aIadapt*(resetV-restV_E)*ones(Ne, 1); 

% --------------------------------------------
% Variables used for learning
sumWee0 = sum(W(1:Ne, 1:Ne), 2); % Initial summed E-to-E weights, for homeostatic normalization

gRiseE = zeros(NetPars.Ncells, 1);
gDecayE = zeros(NetPars.Ncells, 1);
gRiseI = zeros(NetPars.Ncells, 1);
gDecayI = zeros(NetPars.Ncells, 1);

% Initial values for STDP
uSTDP = resetV*zeros(Ne, 1);
vSTDP = resetV*zeros(Ne, 1);
xSTDP = zeros(Ne,1);
yISTDP = zeros(Ncells,1); % for iSTDP

% IdxW:
% 1st col: index of non-zero synaptic weight in array W
% 2nd col: index of presynaptic neuron of corresponding non-zero synaptic connection
% 3nd col: index of postynaptic neuron of corresponding non-zero synaptic connection

% W(i,j) synaptic connection from presynaptic neuron j to postsynaptic neuron i
IdxW = find(W);
[rowSub, colSub] = ind2sub(size(W), IdxW);
IdxW = sparse(rowSub, colSub, IdxW, Ncells, Ncells);

% IdxW = find(W);
% [IdxPost, IdxPre] = ind2sub(size(W), IdxW);
% IdxW = [IdxW, IdxPre, IdxPost];
% 
% IdxW(IdxW(:,3)>Ne, :) = []; % Delete the connections to I neurons, because their connections are fixed
% 
% IdxWee = IdxW( (IdxW(:,2)<=Ne)&(IdxW(:,3)<=Ne), :);
% IdxWei = IdxW( (IdxW(:,2)>Ne)&(IdxW(:,3)<Ne), :);
% clear IdxPre IdxPost

%% Iteration
for iter = 1: nStep
    if mod(iter, nStep/1e2) == 1
       fprintf('%d\n', round(1e2*iter/nStep));
    end        
    %% Calculate synaptic inputs of every neurons
    bSpkFF = (rateExt > rand(NetPars.Ncells, 1)); % feedforward input (binary)
    bSpkFFAarray(:, iter) = bSpkFF;

    % Synaptic input (% External input + recurrent input)
    InFF = zeros(Ncells, 1);
    InFF(bSpkFF) = Jff(bSpkFF) .* bSpkFF(bSpkFF);
    
    InSynE = sum(W(:, bSpkNet(1:Ne)), 2); % Excitatory inputs to all neurons, [N, 1]
    InSynI = sum(W(:, [false(Ne,1); bSpkNet(Ne+1:end)]), 2); % Inhibitory inputs to all neurons, [N, 1]
    %     InSynE = W(:, bSpkNet(1:Ne)) * bSpkNet(bSpkNet(1:Ne));
    %     InSynI = W(:, [false(Ne,1); bSpkNet(Ne+1:end)]) ...
    %         * bSpkNet([false(Ne,1); bSpkNet(Ne+1:end)]);
    
    gRiseE  = (1-dt/tauRiseE)  * gRiseE  + InSynE + InFF;
    gDecayE = (1-dt/tauDecayE) * gDecayE + InSynE + InFF;
    gRiseI  = (1-dt/tauRiseI)  * gRiseI  + InSynI;
    gDecayI = (1-dt/tauDecayI) * gDecayI + InSynI;
    
    gE = (gDecayE - gRiseE) / (tauDecayE - tauRiseE);
    gI = (gDecayI - gRiseI) / (tauDecayI - tauRiseI);
    
    ISyn = gE.*(revV_E - V) + gI.*(revV_I - V);     
    
    % Adaptation current (only for E neurons)
    Iadapt = Iadapt + (aIadapt*(V(1:Ne)-restV_E)-Iadapt)*dt/tauIadapt;
    IA(:, iter) = Iadapt;
    % Dynamic firing threshold for E neurons
    vth = vth + (thetaV - vth)*dt/tauTheta;
    
    %% Updating membrane potential
    % ------------------------
    % Subthreshold activities
    
    % Get the index of neurons which are not in refractory period
    Idx = (iter*dt > (tLastSpk + tauRefrac)); % not in refractory period
    IdxE = Idx(1:Ne);
    IdxI = [false(Ne, 1); Idx(Ne+1: end)];
    
    % Excitatory neurons
    dVe = (restV_E - V(IdxE) + deltaTheta*exp((V(IdxE)-vth(IdxE))/deltaTheta))/tauE ...
        + (ISyn(IdxE) - Iadapt(IdxE))/C;
    V(IdxE) = V(IdxE) + dVe * dt;
    bSpkNet(1:Ne) = (V(1:Ne)>Vpeak);
    
    % Inhibitory neuron
    dVi = (restV_I - V(IdxI))/tauI + ISyn(IdxI)/C;
    V(IdxI) = V(IdxI) + dVi * dt;
    bSpkNet(Ne+1:end) = (V(Ne+1:end)> thetaV);
    
    V(bSpkNet) = resetV; % Reset membrane potential if it fires a spike
    % --------------------------
    % During refractory period
    V(~Idx) = resetV;
    if nargout > 2
        VArray(:,iter+1) = V;
    end
    
    % ---------------------------------------    
    vth(bSpkNet(1:Ne)) = thetaV + ATheta;
    Iadapt(bSpkNet(1:Ne)) = Iadapt(bSpkNet(1:Ne)) + bIadapt;
    
    % --------------------
    % Record spike timing
    tLastSpk(bSpkNet) = iter * dt;   
    if sum(bSpkNet)>0
        tSpk(:, nSpkAll:nSpkAll+sum(bSpkNet)-1) = [find(bSpkNet)'; tLastSpk(bSpkNet)'];
        nSpkAll = nSpkAll+sum(bSpkNet);
    end

  %% Updating the connection weight according to neuronal activities
    % ------------------------------------------
    % Voltage based STDP (modify E-to-E connections)
    % Filtering membrane potential and spike trains
    uSTDP = uSTDP + (V(1:Ne) - uSTDP)*dt/tauu;
    vSTDP = vSTDP + (V(1:Ne) - vSTDP)*dt/tauv;
    xSTDP = xSTDP * (1-dt/taux);
    xSTDP(bSpkNet(1:Ne)) = xSTDP(bSpkNet(1:Ne)) + 1./taux; % When a spike is emitted
    
    % Updating E-to-E connections according to STDP
    IdxWLTD = IdxW(uSTDP>thetaLTD, bSpkNet(1:Ne));
    IdxWLTD = IdxWLTD(IdxWLTD>0);
    if ~isempty(IdxWLTD)
        [IdxPostLTD, ~] = ind2sub(Ncells*ones(1,2), IdxWLTD);
        W(IdxWLTD) = W(IdxWLTD) - aLTD*dt*(uSTDP(IdxPostLTD) - thetaLTD);
        W(IdxWLTD) = max(W(IdxWLTD), JeeMin);
    end
    
    IdxWLTP = IdxW((V(1:Ne)>thetaLTP)&(vSTDP>thetaLTD), 1:Ne);
    IdxWLTP = IdxWLTP(IdxWLTP>0);
    IdxWLTP = IdxWLTP(:);
    if ~isempty(IdxWLTP)
        [IdxPostLTP, IdxPreLTP] = ind2sub(Ncells*ones(1,2), IdxWLTP);
        W(IdxWLTP) = W(IdxWLTP) + ...
            aLTP*dt* xSTDP(IdxPreLTP).*(V(IdxPostLTP)-thetaLTP).*(vSTDP(IdxPostLTP)-thetaLTD);
        W(IdxWLTP) = max(W(IdxWLTP), JeeMax);
    end
        
    % -----------------------------------------------------
    % Homeostatic mechanism (normalize E-to-E connections)
    if mod(iter, dtnormalize/dt) == 0
        dWee = (sum(W(1:Ne, 1:Ne), 2)-sumWee0)./sum(W(1:Ne, 1:Ne)>0, 2);
        W(1:Ne, 1:Ne) = bsxfun(@minus, W(1:Ne, 1:Ne), dWee);
    end
        
    % -------------------------------------------
    % Inhibitory STDP (modify I-to-E connections)
    yISTDP = yISTDP * (1-dt/tauy);
    yISTDP(bSpkNet) = yISTDP(bSpkNet) + 1; % When a neuron fires a spike
    
    % When a presynaptic I neuron fires a spike
    IdxWEI1 = IdxW(:, [false(Ne, 1); bSpkNet(Ne+1:end)]);
    IdxWEI1 = IdxWEI1(IdxWEI1>0);
    if ~isempty(IdxWEI1)
        IdxWEI1 = IdxWEI1(:);
        [IdxPostE, ~] = ind2sub(Ncells*ones(1,2), IdxWEI1);
        W(IdxWEI1) = W(IdxWEI1) + eta*(yISTDP(IdxPostE)-2*r0*tauy);
    end
    
    % When a postesynaptic E neuron fires a spike
    IdxWEI2 = IdxW(bSpkNet(1:Ne), :);
    IdxWEI2 = IdxWEI2(IdxWEI2>0);
    if ~isempty(IdxWEI2)
        IdxWEI2 = IdxWEI2(:);
        [~, IdxPreI] = ind2sub(Ncells*ones(1,2), IdxWEI2);
        W(IdxWEI2) = W(IdxWEI2) + eta*yISTDP(IdxPreI);
    end
    
    IdxWEI = union(IdxWEI1, IdxWEI2);
    W(IdxWEI) = max(IdxWEI, JeiMin);
    W(IdxWEI) = min(IdxWEI, JeiMax);

end