# network_inference_sandbox
Sandbox for toy examples of spiking networks: generating data, and playing with reconstruction algorithms

# Contents

mlsqrnoisy3.m

Matlab code originally from my SVN directory /Users/Shared/projects/neuro/code/morrislecar/matlab
Probably dating back to 2008 or earlier.  Implements noisy Morris-Lecar with specified coupling.
Defaults to three cells but can be more complicated.  

spiketimes.m

Finds the time(s) at which a voltage trace crosses a fixed voltage threshold
vthresh, up to linear accuracy in the time increments, using linear
interpolation.

spiketimeplot.m 

Plots a cell array of spike times as a raster plot.

spiketimesimilarity.m

Calculates a spike time similarity matrix at a given time resolution. 

spiketrains_smoothed.m

Smooths spike trains on different timescales with exponential kernels.

