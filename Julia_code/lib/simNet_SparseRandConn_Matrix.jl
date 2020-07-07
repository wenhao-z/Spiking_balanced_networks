#this file is part of alk_cluster_2012
#Copyright (C) 2014 Ashok Litwin-Kumar
#see README for more information

function simNet_SparseRandConn_matrix(parsNet::Dict{String, Any})
	println("setting up parameters")

	# Ncells::Int64 	= parsNet["Ncells"]
	# Ne::Int64 		= parsNet["Ne"]
	# Ni::Int64 		= parsNet["Ni"]
	# T::Int64		= parsNet["T"]

	Ncells 	= parsNet["Ncells"]
	Ne 		= parsNet["Ne"]
	Ni 		= parsNet["Ni"]
	T		= parsNet["T"]

	#########################################################
	# Intrinsic parameters of single neurons
	#########################################################
	taue = 15 #membrane time constant for exc. neurons (ms)
	taui = 10

	vre = 0. #reset voltage
	threshe = 1 #threshold for exc. neurons
	threshi = 1

	dt = .1 #simulation timestep (ms)
	refrac = 5 #refractory period (ms)

	#synaptic time constants (ms)
	tauerise = 1
	tauedecay = 3
	tauirise = 1
	tauidecay = 2

	#########################################################
	# Parameters of connections
	#########################################################
	# square root of the average number of E->E connections per neuron
	sqrtK = sqrt(parsNet["pee"] * parsNet["Ne"])
	rho = parsNet["Ne"] / (2*parsNet["Width"])

	# Load the connection parameters within and without population
	# jie::Float64 = parsNet["jie_scale"]/(taui*sqrtK)
	# jei::Float64 = parsNet["jei_scale"]/(taue*sqrtK)
	# jii::Float64 = parsNet["jii_scale"]/(taui*sqrtK)
	# jee::Float64 = parsNet["jee_scale"]/(taue*sqrtK)

	jie = parsNet["jie_scale"]/(taui*sqrtK)
	jei = parsNet["jei_scale"]/(taue*sqrtK)
	jii = parsNet["jii_scale"]/(taui*sqrtK)
	jee = parsNet["jee_scale"]/(taue*sqrtK)

	#########################################################
	# Set the connection matrix
	#########################################################
	TuneWidth  = parsNet["TuneWidth"]
	Width = parsNet["Width"]

	Ker = parsNet["PrefStim"] .- parsNet["PrefStim"][1]
	Ker = angle.(exp.(Ker*pi*1im/Width))* Width/pi
	Ker = exp.(-Ker.^2/(2*TuneWidth^2))/(sqrt(2*pi)*TuneWidth)

	W = zeros(parsNet["Ne"], parsNet["Ne"])
	for iter = 1: parsNet["Ne"]
		W[iter,:] = circshift(Ker, iter)
	end
	ratio_jeestruct = parsNet["ratio_jeestruct"]
	W = jee*(ratio_jeestruct*2*Width* W .+ 1 .- ratio_jeestruct)

	weights = zeros(Ncells,Ncells)
	#random connections
	weights[1:Ne,1:Ne] = W .* (rand(Ne,Ne) .< parsNet["pee"])
	weights[1:Ne,(1+Ne):Ncells] = jei*(rand(Ne,Ni) .< parsNet["pei"])
	weights[(1+Ne):Ncells,1:Ne] = jie*(rand(Ni,Ne) .< parsNet["pie"])
	weights[(1+Ne):Ncells,(1+Ne):Ncells] = jii*(rand(Ni,Ni) .< parsNet["pii"])

	# No self connections
	for ci = 1:Ncells
		weights[ci,ci] = 0
	end
	weights_fromE = weights[:, 1:Ne]
	weights_fromI = weights[:, (Ne+1):Ncells]

	#########################################################
	# Set the stimulus/input parameters
	#########################################################
	# Nstim::Int64  = parsNet["Nstim"]
	# stimstr::Float64 = parsNet["stimstr_scale"]/taue
	# stimstart::Int64 = parsNet["stimstart"]
	# stimend::Int64 = parsNet["stimend"]

	Nstim  = parsNet["Nstim"]
	stimstr = parsNet["stimstr_scale"]/taue
	stimstart = parsNet["stimstart"]
	stimend = parsNet["stimend"]

	#constant bias to each neuron type
	muemin = 1.1
	muemax = 1.2
	muimin = 1
	muimax = 1.05

	#########################################################
	# Initialize arrays for network simulation
	#########################################################
	mu = zeros(Ncells)
	mu[1:Ne] = (muemax-muemin)*rand(Ne) .+ muemin
	mu[(Ne+1):Ncells] = (muimax-muimin)*rand(Ni) .+ muimin

	thresh = zeros(Ncells)
	thresh[1:Ne] .= threshe
	thresh[(1+Ne):Ncells] .= threshi

	tau = zeros(Ncells)
	tau[1:Ne] .= taue
	tau[(1+Ne):Ncells] .= taui

	maxTimes = round(Int,parsNet["maxrate"]*T/1000)
	tSpk = zeros(Ncells,maxTimes)
	# tSpk = zeros(2, Ncells*maxTimes)
	nSpk = zeros(Int,Ncells)

	forwardInputsE = zeros(Ncells) #summed weight of incoming E spikes
	forwardInputsI = zeros(Ncells)
	forwardInputsEPrev = zeros(Ncells) #as above, for previous timestep
	forwardInputsIPrev = zeros(Ncells)

	xerise = zeros(Ncells) #auxiliary variables for E/I currents (difference of exponentials)
	xedecay = zeros(Ncells)
	xirise = zeros(Ncells)
	xidecay = zeros(Ncells)
	synInput = zeros(Ncells)

	v = rand(Ncells) #membrane voltage

	lastSpike = -100*ones(Ncells) #time of last spike

	Nsteps = round(Int,T/dt)

	###############################################################
	# Simulation loop
	###############################################################
	println("starting simulation")
	for ti = 1:Nsteps
		if mod(ti,Nsteps/100) == 1  #print percent complete
			@printf("\r%d%%",round(Int,100*ti/Nsteps))
		end
		t = dt*ti
		forwardInputsE[:] .= 0.
		forwardInputsI[:] .= 0.

		# global xerise 	+= -dt*xerise/tauerise 	 + forwardInputsEPrev
		# global xedecay += -dt*xedecay/tauedecay + forwardInputsEPrev
		# global xirise 	+= -dt*xirise/tauirise 	 + forwardInputsIPrev
		# global xidecay += -dt*xidecay/tauidecay + forwardInputsIPrev
		#
		# synInput = (xedecay - xerise)/(tauedecay - tauerise) +
		# 	(xidecay - xirise)/(tauidecay - tauirise)
		# if (t > stimstart) && (t < stimend)
		# 	synInput[1:Nstim] += stimstr;
		# end
		#
		# # idxCell = (t .> (lastSpike .+ refrac))
		# global idxCell = lastSpike.< (t-refrac)
		# v[idxCell] += dt*((1 ./tau[idxCell]).*(mu[idxCell]-v[idxCell]) .+
		# 	synInput[idxCell])
		#
		# idxSpkCell = (v .> thresh)
		# v[idxSpkCell] .= vre
		# lastSpike = t
		# nSpk[idxSpkCell] .+= 1
		#
		# for iter = findall(idxSpkCell)
		# 	tSpk[iter, nSpk[iter]] = t
		# end
		#
		# # Recurrent input
		# forwardInputsE = sum(weights_fromE[:, idxSpkCell[1:Ne]], dims=2)
		# forwardInputsI = sum(weights_fromI[:, idxSpkCell[(Ne+1):Ncells]], dims=2)

		# Recurrent inputs
		xerise  += -dt.*xerise./tauerise   .+ forwardInputsEPrev
		xedecay += -dt.*xedecay./tauedecay .+ forwardInputsEPrev
		xirise  += -dt.*xirise./tauirise   .+ forwardInputsIPrev
		xidecay += -dt.*xidecay./tauidecay .+ forwardInputsIPrev

		synInput = (xedecay - xerise)./(tauedecay - tauerise) +
			(xidecay - xirise)./(tauidecay - tauirise)
		# if (ci < Nstim) && (t > stimstart) && (t < stimend)
		# 	synInput += stimstr;
		# end

		for ci = 1:Ncells
			if t > (lastSpike[ci] + refrac)  #not in refractory period
				v[ci] += dt*((1/tau[ci])*(mu[ci]-v[ci]) + synInput[ci])

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

	# tSpk = tSpk[:,1:maximum(nSpk)]
	# tSpk = tSpk[:,1:min(size(tSpk,2), maximum(nSpk))]

	return tSpk,nSpk, weights
end
