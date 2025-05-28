function S=spiketimesimilarity(spiketime_array,sigma)

%function S=spiketimesimilarity(spiketime_array,sigma)
%
% Calculate the spike time similarity matrix at a given time resolution
% sigma. 
%
% Inputs:
%
% spiketime_array is a cell array with nv spike time records (increasing
% real numbers).
%
% sigma (optional, default = 25 (msec)) represents the time scale of
% comparison.
%
% Outputs:
%
% S is an nv x nv symmetric matrix with the normalized spike time
% similarity between the ith and jth spike trains in Sij.  We assume each
% cell contains a row vector of spike times.
%
% Let g(u)=exp(-u^2/(2 sigma^2))/sqrt(2 pi sigma^2) be the standard
% Guassian with standard deviation sigma.  
%
% Smoothing each spike train by convolution with g(u), we calculate the
% inner product between two a pair of smoothed spike trains analytically,
% and normalize so that entries of S lie between +-1.  
%
% In this version, we do not subtract off the mean firing rate.  
%
% Peter Thomas 2025-05-26 (CWRU), based on the ideas in Fellous, J. M.,
% Tiesinga, P. H., Thomas, P. J., & Sejnowski, T. J. (2004). Discovering
% spike patterns in neuronal responses. Journal of Neuroscience, 24(12),
% 2989-3001.   

if nargin < 2
    sigma=25; % msec
end
nv=length(spiketime_array); % number of spike traces to be compared

% First, find the self-overlap of each spike train, to normalize.
Z=nan(1,nv);
for i=1:nv
    t_i=spiketime_array{i};
    n_i=length(t_i);
    Z(i)=sum(sum(exp(-(ones(n_i,1)*t_i - t_i'*ones(1,n_i)).^2/(4*sigma^2))));
end

% Next, find the overlap between pairs
S=ones(nv);
for i=1:(nv-1)
    t_i=spiketime_array{i};
    n_i=length(t_i);
    for j=(i+1):nv
        t_j=spiketime_array{j};
        n_j=length(t_j);
        S(i,j)=sum(sum(exp(-(ones(n_j,1)*t_i - t_j'*ones(1,n_i)).^2/(4*sigma^2))))/sqrt(Z(i)*Z(j));
        S(j,i)=S(i,j);
    end
end

%% Make a plot of S, if desired
plotting=1;
if plotting
    zvec=zeros(nv,1);
    figure
    pcolor([S,zvec;zvec',0])
    colorbar
    shading flat
    colormap hot
    title('Spike Train Similarity')
    set(gca,'FontSize',20)
end
