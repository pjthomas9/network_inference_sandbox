function spiketime_array=spiketimes(t,v,vthresh)

%function spiketime_array=spiketimes(t,v,vthresh)
%
% For one or several continuous voltage traces sampled at times t, find the
% time(s) at which each trace crosses a fixed voltage threshold vthresh, up
% to linear accuracy in the time increments, using linear interpolation.
%
% Inputs: 
% 
% (1 x nt) vector t
% (nv x nt) array v
% (1 x 1) threshold vthresh
%
% Outputs:
%
% nv component cell array, each cell containing a vector of the spike times
% for a specific voltage trace
%
% PJT 2025-05-26 CWRU

nt=length(t);
[nv,nt2]=size(v);
if nt~=nt2
    warning("Warning! Input v should have number of columns equal to the length of t!")
end
spiketime_array=cell(1,nv); % create cell array to store output

%% Loop through voltage traces
for i=1:nv
    spikeflag=(v(i,1:end-1)<vthresh)&(v(i,2:end)>=vthresh);
    spikeflag1=logical([spikeflag,0]);
    spikeflag2=logical([0,spikeflag]);
    t1=t(spikeflag1);
    v1=v(i,spikeflag1);
    t2=t(spikeflag2);
    v2=v(i,spikeflag2);
    spiketime_array{i}=t1+(t2-t1).*((vthresh-v1)./(v2-v1)); % linear interpolation
end
