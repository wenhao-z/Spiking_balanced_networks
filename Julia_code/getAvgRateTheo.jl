function getAvgRateTheo(parsNet::Dict{Symbol, Any})
    # Calculate the mean firing rate of neurons in spiking neural networks
    # using the method presented in Brunel et al., Neural Computation 1999
    # The method assumes neurons are SPARSELY RANDOM connected.

    # Wen-Hao Zhang, Jan-19, 2018
    # Math Department, University of Pittsburgh

    

    stdSynInput =

    yTheta = threshe - (muemin + muemax)/2
    yRest = threshe - (muemin + muemax)/2
    yTheta = yTheta/stdSynInput
    yRest = yRest/stdSynInput


    muemin::Float64 = 1.1
    muemax::Float64 = 1.2
    muimin::Float64 = 1
    muimax::Float64 = 1.05


end
