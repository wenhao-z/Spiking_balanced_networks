# A demo of network
# Wen-Hao Zhang
# University of Pittsburgh
# Jan-12, 2018

using PyPlot, Printf, LinearAlgebra
# @pyimport numpy as np

# bPlot = false
bSave = false
bPlot = true

# Load default parameters
include("parsNetMdl_default.jl")
include("lib/simNet_SparseRandConn.jl")
# include("lib/simNet_SparseRandConn_matrix.jl")

# Specified parameters
parsNet["T"] = Int(5e3);
parsNet["stimstart"] = 0;
parsNet["stimend"] = 0;

# Perform network simulation
@time tSpk,nSpk,weights = simNet_SparseRandConn(parsNet)
# @time tSpk,nSpk,weights = simNet_SparseRandConn_matrix(parsNet)
# @time tSpk, nSpk, weights = simNetMdl(parsNet)

# include("statSpk.jl")

# Plot figures
if bPlot
	using PyPlot
	println("creating plot")
	PyPlot.clf()
	figure(figsize=(4,4))
	for ci = 1:parsNet["Ne"]
		vals = tSpk[ci,1:nSpk[ci]]
		y = ci*ones(length(vals))
		scatter(vals,y,s=.3,c="k",marker="o",linewidths=0)
	end

	xlim(0,parsNet["T"])
	ylim(0,parsNet["Ne"])
	ylabel("Neuron")
	xlabel("Time")
	tight_layout()
	gcf()

	# savefig("output.png",dpi=150)
end

# Save (HDF5)
if bSave
	using HDF5, Dates
	if ~isdir("Data")
		mkdir("Data")
	end
	timeStr = Dates.format(now(), "yymmdd_HHMM")
	fileName = string("Data/","NetRes_", timeStr, ".h5")
	h5open(fileName, "w") do file
	    write(file, "tSpk", tSpk)
	    write(file, "nSpk", nSpk)
		write(file, "weights", weights)
	end
	h5writeattr(fileName, "tSpk", parsNet)
end
