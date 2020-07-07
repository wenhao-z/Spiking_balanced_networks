#this file is part of alk_cluster_2012
#Copyright (C) 2014 Ashok Litwin-Kumar
#see README for more information

function simNetMdl(parsNet::Dict{String, Real})
	println("setting up parameters")
	# Single neuron parametres
	taue::Int64 	= 15 #membrane time constant for exc. neurons (ms)
	taui::Int64 	= 10
	vre::Float64 	= 0. #reset voltage
	threshe::Int64 	= 1 #threshold for exc. neurons
	threshi::Int64 	= 1
	dt::Float64 	= .1 #simulation timestep (ms)
	refrac::Int64 	= 5 #refractory period (ms)

	#synaptic time constants (ms)
	tauerise::Int64  = 1
	tauedecay::Int64 = 3
	tauirise::Int64  = 1
	tauidecay::Int64 = 2

	#constant bias to each neuron type
	muemin::Float64 = 1.1
	muemax::Float64 = 1.2
	muimin::Float64 = 1
	muimax::Float64 = 1.05

	# Unpack parameters from parsNet, and calculat dependent parameters
	Ncells::Int64 	= parsNet["Ncells"]
	Ne::Int64 		= parsNet["Ne"]
	Ni::Int64 		= parsNet["Ni"]
	Nepop::Int64 	= parsNet["Nepop"]
	T::Int64 		= parsNet["T"]

	# Stimulation
	Nstim::Int64 		= parsNet["Nstim"]
	stimstr::Float64 	= parsNet["stimstr_ratio"]/taue
	stimstart::Int64 	= parsNet["stimstart"]
	stimend::Int64 		= parsNet["stimend"]
	maxrate::Int64 		= parsNet["maxrate"]

	# Connection probability and weight
	sqrtK::Float64 	= sqrt(parsNet["K"])
	jie::Float64 	= parsNet["jie_ratio"]/(taui*sqrtK)
	jei::Float64 	= parsNet["jei_ratio"]/(taue*sqrtK)
	jii::Float64 	= parsNet["jii_ratio"]/(taui*sqrtK)
	jeeout::Float64 = parsNet["jeeout_ratio"]/(taue*sqrtK)
	jeein::Float64 	= parsNet["ratioejee"]*jeeout

	# pei::Float64 	= parsNet[:pei]
	# pie::Float64 	= parsNet[:pie]
	# pii::Float64 	= parsNet[:pii]
	peeout::Float64 = parsNet["K"]/(Nepop*(parsNet["ratiopee"]-1) + Ne)
	peein::Float64 	= parsNet["ratiopee"]*peeout

	#################################################################
	# Initialize variables
	mu = zeros(Ncells)
	mu[1:Ne] = (muemax-muemin)*rand(Ne) + muemin
	mu[(Ne+1):Ncells] = (muimax-muimin)*rand(Ni) + muimin

	thresh = zeros(Ncells)
	thresh[1:Ne] = threshe
	thresh[(1+Ne):Ncells] = threshi

	tau = zeros(Ncells)
	tau[1:Ne] = taue
	tau[(1+Ne):Ncells] = taui

	Npop = round(Int,Ne/Nepop)

	# Initialize weight matrix (random connections)
	weights = zeros(Ncells,Ncells)
	weights[1:Ne,1:Ne] = jeeout*(rand(Ne,Ne) .< peeout)
	weights[1:Ne,(1+Ne):Ncells] = jei*(rand(Ne,Ni) .< parsNet["pei"])
	weights[(1+Ne):Ncells,1:Ne] = jie*(rand(Ni,Ne) .< parsNet["pie"])
	weights[(1+Ne):Ncells,(1+Ne):Ncells] = jii*(rand(Ni,Ni) .< parsNet["pii"])

	#connections within cluster
	for iter_pop = 1:Npop
		ipopstart = 1 + Nepop*(iter_pop-1)
		ipopend = iter_pop*Nepop

		weights[ipopstart:ipopend,ipopstart:ipopend] = jeein*(rand(Nepop,Nepop) .< peein)
	end

	for ci = 1:Ncells
		weights[ci,ci] = 0
	end

	maxTimes 	= round(Int,maxrate*T/1000)
	tSpk 		= zeros(Ncells,maxTimes)
	nSpk 		= zeros(Int,Ncells)

	forwardInputsE 		= zeros(Ncells) #summed weight of incoming E spikes
	forwardInputsI 		= zeros(Ncells)
	forwardInputsEPrev 	= zeros(Ncells) #as above, for previous timestep
	forwardInputsIPrev 	= zeros(Ncells)

	xerise 	= zeros(Ncells) #auxiliary variables for E/I currents (difference of exponentials)
	xedecay = zeros(Ncells)
	xirise 	= zeros(Ncells)
	xidecay = zeros(Ncells)

	v = rand(Ncells) #membrane voltage

	lastSpike = -100*ones(Ncells) #time of last spike

	Nsteps = round(Int,T/dt)

	println("starting simulation")

	#begin main simulation loop
	for ti = 1:Nsteps
		if mod(ti,Nsteps/10) == 1  #print percent complete
			@printf("\r%d%%",round(Int,100*ti/Nsteps))
		end
		t = dt*ti
		forwardInputsE[:] = 0
		forwardInputsI[:] = 0

		for ci = 1:Ncells
			xerise[ci] 	+= -dt*xerise[ci]/tauerise + forwardInputsEPrev[ci]
			xedecay[ci] += -dt*xedecay[ci]/tauedecay + forwardInputsEPrev[ci]
			xirise[ci] 	+= -dt*xirise[ci]/tauirise + forwardInputsIPrev[ci]
			xidecay[ci] += -dt*xidecay[ci]/tauidecay + forwardInputsIPrev[ci]

			synInput = (xedecay[ci] - xerise[ci])/(tauedecay - tauerise) +
			(xidecay[ci] - xirise[ci])/(tauidecay - tauirise);

			if (ci < Nstim) && (t > stimstart) && (t < stimend)
				synInput += stimstr;
			end

			if t > (lastSpike[ci] + refrac)  #not in refractory period
				v[ci] += dt*((1/tau[ci])*(mu[ci]-v[ci]) + synInput)

				if v[ci] > thresh[ci]  #spike occurred
					v[ci] = vre
					lastSpike[ci] = t
					nSpk[ci] = nSpk[ci]+1
					if nSpk[ci] <= maxTimes
						tSpk[ci,nSpk[ci]] = t
					end

					for j = 1:Ncells
						if weights[j,ci] > 0  #E synapse
							forwardInputsE[j] += weights[j,ci]
						elseif weights[j,ci] < 0  #I synapse
							forwardInputsI[j] += weights[j,ci]
						end
					end #end loop over synaptic projections
				end #end if(spike occurred)
			end #end if(not refractory)
		end #end loop over neurons

		forwardInputsEPrev = copy(forwardInputsE)
		forwardInputsIPrev = copy(forwardInputsI)
	end #end loop over time
	@printf("\r")

	tSpk = tSpk[:,1:maximum(nSpk)]

	return tSpk,nSpk,weights
	# return tSpk,nSpk,Ne,Ncells,T
end
