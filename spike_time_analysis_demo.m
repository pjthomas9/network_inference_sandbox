% Script to illustrate use of the tools in this directory.

% In this illustration, we have already run the following code:
%
% c=ones(3)-eye(3); % this is a self-connected excitatory block
% cmat=[c,-ones(3),zeros(3);zeros(3),c,-ones(3);-ones(3),zeros(3),c];
% (Remember to set initial conditions with init = '3-3-3-antisync' in mlsqrnoisy3.m;)
% [t,y,inputs]=mlsqrnoisy3(0.2,2e4-1,cmat,1e-2/3,1e-2/6,1e3,1);
%
% This code simulates nine coupled stochastic Morris Lecar neurons with
% groups of three coupled in mutually excitatory blocks, and the ith block
% inhibiting the i+1st block to produce a three-phase pattern.
%
% The output was stored in mlsqrn3_20250525T214322.mat, which we load
% first:

load('mlsqrn3_20250525T214322.mat') % loads inputs, t, and y.

% Extract the voltage 
nv=9; % number of voltage traces
v=y(1:nv,:);

% Extract the spike times
vthresh=4.3; % empirically chosen
spiketime_array=spiketimes(t,v,vthresh);

% Plot the spike times

tmin=0;
tmax=2e4;
[h,idx]=spiketimeplot(spiketime_array,tmin,tmax);

% Calculate a spike time similarity measure at a specified resolution
sigma=25; % msec
% Shuffle spike train labels
spiketime_array_shuffled=cell(size(spiketime_array));
for i=1:nv
    spiketime_array_shuffled{i}=spiketime_array{idx(i)};
end
S_shuffled=spiketimesimilarity(spiketime_array_shuffled,sigma);
