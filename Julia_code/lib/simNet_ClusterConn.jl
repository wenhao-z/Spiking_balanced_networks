#this file is part of alk_cluster_2012
#Copyright (C) 2014 Ashok Litwin-Kumar
#see README for more information

function simNet_SparseRandConn(parsNet::Dict{String, Real})
	println("setting up parameters")

	Ncells::Int64 	= parsNet["Ncells"]
	Ne::Int64 		= parsNet["Ne"]
	Ni::Int64 		= parsNet["Ni"]
	Nepop::Int64 	= parsNet["Nepop"]
	T::Int64		= parsNet["T"]

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
	# Load connection probabilities
	pee::Float64  = parsNet["pee"]
	pei::Float64  = parsNet["pei"]
	pie::Float64  = parsNet["pie"]
	pii::Float64  = parsNet["pii"]

	K = pee * Ne #average number of E->E connections per neuron
	sqrtK = sqrt(K)

	ratiopee::Float64 = parsNet["ratiopee"]
	peeout = K/(Nepop*(ratiopee-1) + Ne)
	peein = ratiopee*peeout

	# Load the connection parameters within and without population
	jie::Float64 	= parsNet["jie_scale"]/(taui*sqrtK)
	jei::Float64 	= parsNet["jei_scale"]/(taue*sqrtK)
	jii::Float64 	= parsNet["jii_scale"]/(taui*sqrtK)
	jeeout::Float64 = parsNet["jeeout_scale"]/(taue*sqrtK)
	jeein::Float64 	= parsNet["ratioejee"]*jeeout

	#########################################################
	# Set the connection matrix
	#########################################################
	Npop = round(Int,Ne/Nepop)

	weights = zeros(Ncells,Ncells)
	#random connections
	weights[1:Ne,1:Ne] = jeeout*(rand(Ne,Ne) .< peeout)
	weights[1:Ne,(1+Ne):Ncells] = jei*(rand(Ne,Ni) .< pei)
	weights[(1+Ne):Ncells,1:Ne] = jie*(rand(Ni,Ne) .< pie)
	weights[(1+Ne):Ncells,(1+Ne):Ncells] = jii*(rand(Ni,Ni) .< pii)

	#connections within cluster
	for IdxNpop = 1:Npop
		ipopstart = 1 + Nepop*(IdxNpop-1)
		ipopend = IdxNpop*Nepop
		weights[ipopstart:ipopend,ipopstart:ipopend] = jeein*(rand(Nepop,Nepop) .< peein)
	end
	# No self connections
	for ci = 1:Ncells
		weights[ci,ci] = 0
	end

	#########################################################
	# Set the stimulus/input parameters
	#########################################################
	Nstim::Int64  = parsNet["Nstim"]
	stimstr::Float64 = parsNet["stimstr_scale"]/taue
	stimstart::Int64 = parsNet["stimstart"]
	stimend::Int64 = parsNet["stimend"]

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
	nSpk = zeros(Int,Ncells)

	forwardInputsE = zeros(Ncells) #summed weight of incoming E spikes
	forwardInputsI = zeros(Ncells)
	forwardInputsEPrev = zeros(Ncells) #as above, for previous timestep
	forwardInputsIPrev = zeros(Ncells)

	xerise = zeros(Ncells) #auxiliary variables for E/I currents (difference of exponentials)
	xedecay = zeros(Ncells)
	xirise = zeros(Ncells)
	xidecay = zeros(Ncells)

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

		for ci = 1:Ncells
			xerise[ci] += -dt*xerise[ci]/tauerise + forwardInputsEPrev[ci]
			xedecay[ci] += -dt*xedecay[ci]/tauedecay + forwardInputsEPrev[ci]
			xirise[ci] += -dt*xirise[ci]/tauirise + forwardInputsIPrev[ci]
			xidecay[ci] += -dt*xidecay[ci]/tauidecay + forwardInputsIPrev[ci]

			synInput = (xedecay[ci] - xerise[ci])/(tauedecay - tauerise) +
				(xidecay[ci] - xirise[ci])/(tauidecay - tauirise)

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

	return tSpk,nSpk, weights
end
