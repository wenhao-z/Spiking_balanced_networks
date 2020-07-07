# Get the statistics of spike trains
# Wen-Hao Zhang
# University of Pittsburgh
# Feb-1, 2018

# Load previous tSpk
using HDF5, PyPlot, StatsBase
fileName = string("Data/","NetRes_180202_2213.h5")
# fileName = string("Data/","NetRes_180202_1641.h5")

tSpk = h5read(fileName, "tSpk")
nSpk = h5read(fileName, "nSpk")
parsNet = h5readattr(fileName, "tSpk")

# Run a long single trial and partition it into several small segments
tSeg = Int(2e3) # 2 seconds
tEdge = collect(0:tSeg:parsNet["T"])

nSpk_trial = zeros(parsNet["Ncells"], length(tEdge)-1)
for iter = 1: parsNet["Ncells"]
    nSpk_trial[iter,:] = fit(Histogram, tSpk[iter,:], tEdge, closed=:left).weights
end

# ------------------------------------------------------------------------------
# Histogram of mean firing rate
rateAvg = nSpk / (parsNet["T"]/1e3)
rEdge = linspace(0, ceil(maximum(rateAvg)), 50)
# rEdge = 0:1:ceil(maximum(rateAvg))
histRate = fit(Histogram, rateAvg, rEdge, closed=:left).weights
subplot(3,2,1)
step((rEdge[1:end-1] + rEdge[2:end])/2, histRate)
xlabel("Mean rate (Hz)")
ylabel("Count")
tight_layout()

# ------------------------------------------------------------------------------
# Fano factor
varRate = var(nSpk_trial[:,2:end]/(tSeg/1e3), 2)
meanRate = mean(nSpk_trial[:,2:end]/(tSeg/1e3), 2)
fanoFactor = varRate./rateAvg
fanoFactor = fanoFactor[~isnan(fanoFactor)]

# fEdge = 0:0.05:5
fEdge = linspace(0, ceil(maximum(fanoFactor)), 50)
histFF = fit(Histogram, vec(fanoFactor), fEdge, closed=:left).weights
subplot(3,2,2)
step((fEdge[1:end-1] + fEdge[2:end])/2, histFF)
xlabel("Fano factor")
ylabel("Count")
# var # variance
# cov # covariance
# cor # Pearson correlation coefficient


# ------------------------------------------------------------------------------
# Correlation coefficients of neuronal firing rate
corSpkMat = cor(nSpk_trial[:, 2:end]')
corSpk = corSpkMat[~isnan(corSpkMat[:])];
cEdge = linspace(floor(minimum(corSpk[:])), ceil(maximum(corSpk[:])), 50)
histCorSpk = fit(Histogram, vec(corSpk), cEdge, closed=:left).weights
subplot(3,2,3)
step((cEdge[1:end-1]+cEdge[2:end])/2, histCorSpk)


corSpk_samepop = zeros(parsNet["Nepop"], parsNet["Nepop"], parsNet["Ne"]/parsNet["Nepop"])
for iter = 1: size(corSpk_samepop)[3]
    corSpk_samepop[:,:,iter] = corSpkMat[(iter-1)*parsNet["Nepop"]+1: iter*parsNet["Nepop"],
                (iter-1)*parsNet["Nepop"]+1: iter*parsNet["Nepop"]]
    corSpk_samepop[find(eye(size(corSpk_samepop)[1]))] = NaN
end
corSpk_samepop = corSpk_samepop[~isnan(corSpk_samepop[:])]
cEdge = linspace(floor(minimum(corSpk[:])), ceil(maximum(corSpk[:])), 50)
histCorSpk_samepop = fit(Histogram, vec(corSpk_samepop), cEdge, closed=:left).weights
subplot(3,2,4)
step((cEdge[1:end-1]+cEdge[2:end])/2, histCorSpk_samepop)

# ------------------------------------------------------------------------------
# Cross-correlation function
tBin = 2
tEdge = collect(0:tBin:parsNet["T"])
tlag = -100:100 # unit: number of bins

bSpk = zeros(parsNet["Ncells"], length(tEdge)-1)
for iter = 1: parsNet["Ncells"]
    bSpk[iter,:] = fit(Histogram, tSpk[iter,:], tEdge, closed=:left).weights
end

# Mean cross correlation function between neurons within the same cluster
meanXcorr_samepop = zeros(length(tlag))
nCount = 1
for iter1 = 1: parsNet["Nepop"]
    for iter2 = 1:iter1-1
        # Binarize tSpk
        xcorrSpk = crosscor(bSpk[iter1,:], bSpk[iter2,:], tlag, demean=false)
        meanXcorr_samepop = (meanXcorr_samepop * nCount + xcorrSpk)/ (nCount + 1)
        nCount += 1
    end
end

meanXcorr_diffpop = zeros(length(tlag))
nCount = 1
for iter1 = 1: parsNet["Nepop"]
    for iter2 = parsNet["Nepop"]+1: 2*parsNet["Nepop"]
        # Binarize tSpk
        xcorrSpk = crosscor(bSpk[iter1,:], bSpk[iter2,:], tlag, demean=false)
        meanXcorr_diffpop = (meanXcorr_diffpop * nCount + xcorrSpk)/ (nCount + 1)
        nCount += 1
    end
end

subplot(3,2,5)
step(tlag*tBin, meanXcorr_samepop)
xlabel("Lag (ms)")
title("Cross-corr (same clusters)")
tight_layout()

subplot(3,2,6)
step(tlag*tBin, meanXcorr_diffpop)
xlabel("Lag (ms)")
ylabel("correlation")
title("Cross-corr (different clusters)")
tight_layout()
