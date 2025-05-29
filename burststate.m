function b=burststate(t,spiketime_array,isithresholds)

%function b=burststate(t,spiketime_array,isithresholds)
%
% Given an array of nv spike trains, classify each point in time as "burst"
% or "not burst".  
%
% Input:
%
% t -- vector of times underlying the original continuous voltage traces
%
% spiketime_array -- a cell array of nv spike trains extracted from nv
% voltage traces, for example using spiketimes.m
%
% isithresholds -- Interspike Interval (ISI) threshold defining the
% difference between an isolated spike and a burst.  Default
% isithresholds=500 msec.  If all cells have the same ISI threshold, then
% isithresholds can be a single number.  Otherwise isithreshold should be a
% vector of nv components (one threshold for each spike train).  msec
%
% PJT 2025-05-28 CWRU

% In sandbox mode, try for a single spike train first:

% for multiple spike trains, expand isithresholds into a vector as needed
if nargin < 3, isithresholds=500; end % msec

st=spiketime_array{1};
isi=diff(st); % interspike intervals

%% Find the burst state at each point in time

burst_state=zeros(size(t));

% Long ISI followed by short ISI marks start of a burst:
burst_init_times=st((isi(1:end-1)>=isithresholds(1))&(isi(2:end)<500));

% Short ISI followed by long ISI marks termination of a burst:
burst_term_times=st((isi(1:end-1)<isithresholds(1))&(isi(2:end)>=500));

nburst=length(burst_init_times);

for i=1:nburst-1
    burst_state=burst_state+(t>=burst_init_times(i)).*(t<burst_term_times(i));
end

keyboard

b=burst_state;

figure,plot(t,v(1,:)),hold on,plot(t,-50+65*burst_state)

% Didn't work. Try this instead:

figure,plot(t,v(1,:)),hold on,plot(st((isi(1:end-1)<500)),4.3,'r+'),plot(st((isi(2:end)<500)),4.3,'r+')
