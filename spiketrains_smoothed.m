function u=spiketrains_smoothed(t,v,vthresh,tau_vec)

% function u=spiketrains_smoothed(t,v,vthresh,tau_vec)
%
% Smooth spike trains on different timescales with exponential kernel.
%
% Inputs:
%
% t -- time vector (in msec) underlying continuous voltage trace
%
% v -- nv x nt array of voltages (in mV), where nv is the number of voltage
% traces (number of rows) and nt is the length of t (number of columns)
%
% tau_vec -- an array of one or several time constants for filtering.
% Default: [4 20 100 500 2500] (msec)
%
% Outputs:
%
% u -- (nv*ntau)) x nt array of filtered voltages, where
% ntau=length(tau_vec). The first ntau rows are v(1,:) convolved with
% exp(-t/tau)*(t>=0) for each value of tau in tau_vec.  

if nargin < 4, tau_vec=[4 20 100 500 2500]; end % Default time constants
if nargin < 3, vthresh=4.3; end % Default voltage threshold (mV)

nt=length(t); % number of time points
nv=size(v,1); % number of voltage traces
ntau=length(tau_vec); % number of time constants for filtering
u=nan(nv*ntau,nt-1);

% Create a bank of filters
filt=zeros(ntau,2*nt-1);
for i=1:ntau
    filt(i,nt:end)=exp(-t/tau_vec(i));
end

% Plot the filters, for illustration.
%figure,for i=1:ntau,subplot(ntau,1,i),plot([-t(end:-1:2),t],filt(i,:)),set(gca,'FontSize',16),ylabel(['\tau=',num2str(tau_vec(i))]),end,xlabel('Time (msec)'),subplot(ntau,1,1),title('Column of Spike Time Filters')

%% Convolve the spike trains with decaying exponentials

for i=1:ntau
    for j=1:nv
        u(i+ntau*(j-1),:)=conv(((v(1,1:(end-1))<vthresh).*(vthresh<=v(1,2:end))),filt(i,:),'same');
    end
end

%% Plot results for one of the spike trains

ntau=length(tau_vec);
figure
subplot(ntau+1,1,1)
plot(t/1000,v(1,:))
hold on
plot(t/1000,vthresh*ones(size(t)))
set(gca,'FontSize',16)
ylabel('V(t)')
title('Spike train convolved with exponentials')
%xlim([0 .6])
grid on
for i=1:ntau
    subplot(ntau+1,1,i+1)
    plot(t(1:end-1)/1000,u(i,:))
    ylabel(['\tau=',num2str(tau_vec(i))])
    set(gca,'FontSize',16)
    %xlim([0 .6])
    grid on
end
xlabel('Time (sec)')
