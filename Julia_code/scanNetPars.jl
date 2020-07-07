# Simulate some simple experiments on spiking neural networks
# Wen-Hao Zhang
# Jan-17, 2018

# Load default parameters
# include("parsNetMdl_default.jl")
include("parsNetMdl_default.jl")

# Some particular parameters
parsNet[:ratioejee]  = 1.
parsNet[:ratiopee]   = 1.
parsNet[:T] = Int(2e3)

parsNet[:Nstim] = Int(4e3)
parsNet[:stimstart] = 100
# parsNet[:stimend] = parsNet[:T]

parsNet[:stimstr_ratio] = .07

################################################################################
# Generate parameter grid
include("lib/genParsGrid.jl")

parsRangeDict = Dict(:stimstr_ratio => collect(0: 0.1:1))

parsGridNet = genParsGrid(parsRangeDict, parsNet)

################################################################################
# Perform network simulation
nSpkE = zeros(size(parsGridNet))
nSpkI = zeros(size(parsGridNet))

include("simNetMdl.jl")

tic()
for iter = 1: length(parsGridNet)

    # Perform network simulation
    # @time ~, nSpk = simNetMdl(parsGridNet[iter])
    ~, nSpk = simNetMdl(parsGridNet[iter])
    nSpkE[iter] = mean(1000*nSpk[1:parsNet[:Ne]]/parsNet[:T])
end
toc()
# Save parameters

# Statistics
# @unpack Ne,Ncells,T = parsNet
# println("mean excitatory firing rate: ",mean(1000*nSpk[1:Ne]/T)," Hz")
# println("mean inhibitory firing rate: ",mean(1000*nSpk[(Ne+1):Ncells]/T)," Hz")
