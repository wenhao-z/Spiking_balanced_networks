#this is how to save compressed connectivity data in the HDF5 format
fid = h5open("data.h5","w")

g = g_create(fid,"data")

g["popmembers","chunk", size(popmembers),"compress",9] = popmembers
g["weights","chunk",size(weights),"compress",9] = weights

close(fid)
