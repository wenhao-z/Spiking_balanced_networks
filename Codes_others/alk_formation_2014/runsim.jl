#this file is part of alk_formation_2014
#Copyright (C) 2014 Ashok Litwin-Kumar
#see README for more information

#uncomment the line below and set doplot=true to plot a raster
using PyPlot
doplot = true
# doplot = false

#if true, loads trained network.  if false, generates a new network.
#loading trained network requires the HDF5 package (to load it, uncomment the line below)
#using HDF5
loadtrained = false

include("sim.jl")

#stimulus matrix is in the following format:
#column 1: index of stimulated population
#columns 2/3: start/stop time of stimulus
#column 4: rate of external drive (in kHz)
stim = zeros(1,4)
stim[1,1] = 1

tic()
if loadtrained
	fid = h5open("trained.h5","r")
	popmembers = read(fid["data"]["popmembers"])
	weights = read(fid["data"]["weights"])
	close(fid)

	times,ns,Ne,Ncells,T = sim(stim,weights,popmembers)
else
	times,ns,popmembers,Ne,Ncells,T,VArray = simnew(stim)
end
toc()

println("mean excitatory firing rate: ",mean(1000*ns[1:Ne]/T)," Hz")
println("mean inhibitory firing rate: ",mean(1000*ns[(Ne+1):Ncells]/T)," Hz")

if doplot
	Npop = size(popmembers,1)
	Nmaxmembers = size(popmembers,2)
	println("creating plot")
	figure(figsize=(4,4))
	xlim(0,T)
	ylim(0,sum(popmembers.>0))
	ylabel("Neuron")
	xlabel("Time")
	tight_layout()

	#plot raster with the order of rows determined by population membership
	rowcount = 0
	for pp = 1:20
		@printf("\rpopulation %d",pp)
		for cc = 1:Nmaxmembers
			if popmembers[pp,cc] < 1
				break
			end
			rowcount+=1
			ind = popmembers[pp,cc]
			vals = times[ind,1:ns[ind]]
			y = rowcount*ones(length(vals))
			scatter(vals,y,s=.3,c="k",marker="o",linewidths=0)
		end
	end
	@printf("\rdone creating plot")
	savefig("output.png",dpi=150)
end
