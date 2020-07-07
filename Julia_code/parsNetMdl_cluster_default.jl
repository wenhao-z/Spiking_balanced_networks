# Set the default network Parameters
# The parameters are packed into a Dictionary

# Wen-Hao Zhang
# Math Department, University of Pittsburgh
# Jan-18, 2018

parsNet = Dict(
"Ncells"  => 5000,
"Ne"      => 4000,
"Ni"      => 1000,
"Nepop" => 80, # the number of neurons in each cluster

#connection probabilities
"pee" => .2,
"pei" => .5,
"pie" => .5,
"pii" => .5,
"ratiopee" => 1, # peein/peeout, the ratio of connetion probability within and between clusters
# "ratiopee" => 2.5, # peein/peeout, the ratio of connetion probability within and between clusters

# Relative connection strength
"jie_scale" => 4.,
"jei_scale" => -16. *1.2,
"jii_scale" => -16.,
"jeeout_scale"   => 10.,
# "ratioejee"      => 1.9, # jeein/jeeout
"ratioejee"      => 1, # jeein/jeeout

"T"              => Int(2e3), #simulation time (ms)
"Nstim"          => 400, # The number of neurons receive feedforward inputs
"stimstr_scale"  => .07, # The relative value of stimulus strength (intensity)
"stimstart"      => 500, # (ms)
"stimend"        => Int(2e3), # (ms)

"maxrate" => 100 #(Hz) maximum average firing rate.
# if the average firing rate across the simulation for any neuron exceeds
# this value, some of that neuron's spikes will not be saved
)
