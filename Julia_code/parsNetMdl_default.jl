# Set the default network Parameters
# The parameters are packed into a Dictionary

# Wen-Hao Zhang
# Math Department, University of Pittsburgh
# Jan-18, 2018

parsNet = Dict{String, Any}(
"Ncells"  => 5000,
"Ne"      => 4000,
"Ni"      => 1000,
"Nepop" => 80, # the number of neurons in each cluster

"Width" => 180, # The width of feature space is 2*width, unit: degree
"TuneWidth" => 40, # 40 degrees

#connection probabilities
"pee" => .2,
"pei" => .5,
"pie" => .5,
"pii" => .5,

# Relative connection strength
"jie_scale" => 4.,
"jei_scale" => -16. *1.2,
"jii_scale" => -16.,
"jee_scale" => 10.,
"ratio_jeestruct" => 0, # the ratio of all structured E-E conns over all E-E conns

"T"              => Int(2e3), #simulation time (ms)
"Nstim"          => 400, # The number of neurons receive feedforward inputs
"stimstr_scale"  => .07, # The relative value of stimulus strength (intensity)
"stimstart"      => 500, # (ms)
"stimend"        => Int(2e3), # (ms)

"maxrate" => 150 #(Hz) maximum average firing rate.
# if the average firing rate across the simulation for any neuron exceeds
# this value, some of that neuron's spikes will not be saved
)


parsNet["PrefStim"] = range(-parsNet["Width"], parsNet["Width"],
    length = parsNet["Ne"]+1)
parsNet["PrefStim"] = parsNet["PrefStim"][2:end]

# TuneWidth  = parsNet["TuneWidth"]
# Width = parsNet["Width"]
# PrefStim = range(-Width, Width, length = parsNet["Ne"]+1)
# PrefStim = PrefStim[2:end]
#
# # a = exp.(PrefStim*pi*im/Width)
#
# Ker = angle.(exp.((PrefStim.- PrefStim[1])*pi*1im/Width))* Width/pi
# Ker = exp.(-Ker.^2/(2*TuneWidth^2))/(sqrt(2*pi)*TuneWidth)
#
# W = zeros(parsNet["Ne"], parsNet["Ne"])
# for iter = 1: parsNet["Ne"]
#     W[iter,:] = circshift(Ker, iter)
# end
