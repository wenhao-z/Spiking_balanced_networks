function NetPars = initConnMat(NetPars)
% Initialize the connection matrix in spiking neural network
% Wen-Hao Zhang, Oct-26, 2016
% wenhaoz1@andrew.cmu.edu
% @Carnegie Mellon University


W = [NetPars.Jee0 * ones(NetPars.Ne), NetPars.Jei0 * ones(NetPars.Ne, NetPars.Ni);
    NetPars.Jie * ones(NetPars.Ni, NetPars.Ne), NetPars.Jii * ones(NetPars.Ni)];

W = W.* (rand(NetPars.Ncells) < NetPars.connProb);

% No self connections
W = W - diag(diag(W));

NetPars.W = sparse(W);