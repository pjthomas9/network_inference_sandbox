function [h,idx]=spiketimeplot(spiketime_array,tmin,tmax)

%function h=spiketimeplot(spiketime_array,tmin,tmax)
%
% Plot a cell array of spike times as a raster plot
%
% Input:
%
% spiketime_array is a cell array with nv entries, (one for each voltage
% trace).  Each entry contains an increasing array of spike times.
% 
% tmin (optional) < tmax (optional) are the max and min times for plotting
%
% Output:
%
% h, a handle to the figure produced
%
% idx, a permutation of the integers 1 to nv.  If viewing the data in the
% original order, idx=1:nv.  As an option, one can generate a randomly
% shuffled display order.  

nv=length(spiketime_array); 

% Uncomment one of these two lines:
idx=1:9; % In this case, don't shuffle the spike trains
%idx=randperm(nv); % In this case, shuffle the spike trains

% find the largest last spike time
tmax_temp=0;
for i=1:nv
    tmax_temp=max(tmax_temp,max(spiketime_array{i}));
end
if nargin < 3
    tmax=1.05*tmax_temp;
end
clear tmax_temp
if nargin < 2
    tmin=0;
end

h=figure;
axis([tmin tmax 0.1 nv+.9])
hold on
for i=1:nv
    t_temp=spiketime_array{idx(i)};
    for j=1:length(t_temp)
        tt=t_temp(j);
        set(line([tt,tt],i+0.4*[-1,1]),'LineWidth',2,'Color','k')
    end
end
xlabel('Time (msec)')
ylabel('Spike Train')
set(gca,'FontSize',20)
title('Spike Time Raster Plot')

