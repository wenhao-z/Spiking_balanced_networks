% Perform the linear stability analysis of network

% Wen-Hao Zhang,
% wenhao.zhang@pitt.edu
% University of Pittsburgh
% Mar 12, 2019

parsNet.T = 2e3; 

parsNet.jxe = 1; % E synaptic weight
parsNet.ratiojie = 4; % The ratio between I synapse over E synapse (absolute value)
parsNet.ratio_jeestruct = 0;
parsNet.FanoFactorInput = 0.01;

parsNet.AmplIff = 0;
parsNet.IBkg_scale = [0.8; 1];

% Calculate dependent parameters
parsNet = getDependentPars(parsNet);

%% Run network simulation
% Optional output arguments

tic
outSet = simSpkNet(parsNet);
toc

fireRate = accumarray(outSet.tSpk(1,:)', 1)/parsNet.T*1e3;
fprintf('Mean rate of E Neurons: %3.2fHz \n', mean(fireRate(1:parsNet.Ne)))
fprintf('Mean rate of I Neurons: %3.2fHz \n', mean(fireRate(parsNet.Ne+1:end)))

%% Linear stability analysis 
addpath(fullfile(Path_RootDir, 'TheoreticalAnalysis'));

% Find the equilibrium state
[rateTheory, Isyn, VarIsyn] = getFiringRateTheory(parsNet);

% Find the derivative of IO function on equilibrium point
dIsyn = 1e-5;
dVarIsyn = 1e-5;

% Perturb the mean of synaptic input
ratePertub = zeros(2,1);
ratePertub(1) = nlFunc_MeanField(Isyn(1)+dIsyn, sqrt(VarIsyn(1)), parsNet.threshe, parsNet.taue, parsNet);
ratePertub(2) = nlFunc_MeanField(Isyn(2)+dIsyn, sqrt(VarIsyn(2)), parsNet.threshi, parsNet.taui, parsNet);

drate_dIsyn = (ratePertub*1e3 - rateTheory)/ dIsyn;

% Perturb the variance of synaptic input
ratePertub = zeros(2,1);
ratePertub(1) = nlFunc_MeanField(Isyn(1), sqrt(VarIsyn(1)+dVarIsyn), parsNet.threshe, parsNet.taue, parsNet);
ratePertub(2) = nlFunc_MeanField(Isyn(2), sqrt(VarIsyn(2)+dVarIsyn), parsNet.threshi, parsNet.taui, parsNet);

drate_dVarIsyn = (ratePertub*1e3 - rateTheory)/ dVarIsyn;

clear ratePertub dIsyn dVarIsyn

% ----------------------------------------------------------------------
% Stability (Jacobian) matrix
weights = genNetConnMat(parsNet);

% dF/dIsyn and dF/dVarIsyn
% The reason to divide by 1e3 in the end is because the tau has unit of ms
gradIsyn = [drate_dIsyn(1) * parsNet.taue * ones(parsNet.Ne,1); ...
    drate_dIsyn(2) * parsNet.taui * ones(parsNet.Ni,1)]/1e3;
gradVarIsyn = [drate_dVarIsyn(1) * parsNet.taue * ones(parsNet.Ne,1); ...
    drate_dVarIsyn(1) * parsNet.taui * ones(parsNet.Ni,1)]/1e3;

JacobMat = gradIsyn .* weights + gradVarIsyn .* weights;

% Eigen-analysis of stability matrix
[eigVec_JacobMat, eigVal_JacobMat] = eig(JacobMat);

eigVal_JacobMat = diag(eigVal_JacobMat);

%% Plot the eigen-spectrum
figure;

subplot(1,2,1)
plot(real(eigVal_JacobMat), imag(eigVal_JacobMat), '.')
hold on
axisLim = 1.5;
plot(ones(1,2), axisLim*[-1,1], '--k')
% axis(axisLim*[-1, 1, -1, 1])
axis square

subplot(1,2,2)
radius_eigVal = abs(eigVal_JacobMat);
[histEigVal, edgeEigVal] = histcounts(radius_eigVal, 30);
stairs((edgeEigVal(1:end-1)+edgeEigVal(2:end))/2, histEigVal)
axis square

