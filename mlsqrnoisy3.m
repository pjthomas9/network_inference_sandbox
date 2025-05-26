function [t,y,inputs]=mlsqrnoisy3(dt0,tmax,cmat,ge,gi,Nm,plotflag,rseed)

%function [t,y,inputs]=mlsqrnoisy3(dt0,tmax,cmat,ge,gi,Nm,plotflag,rseed)
%
% An arbitrary number (typically 3) of coupled
% square-wave bursters made of Morris-Lecar cell
% from http://www.math.pitt.edu/~bard/xpp/odes/mlsqr.ode
% adapted for matlab.
%
% Coupling added in the form suggested in Fall et al, Computational Cell
% Biology, chapter 6.
%
% dt0 = baseline timestep for variable timestep (constant d(arc length))
% tmax = maximum time
% cmat = coupling matrix.  cmat(1,2) is coupling from cell 2 to 1, 1 for
% excitatory and -1 for inhibitory and 0 for no coupling.  
%
% Format for multiple cells: v, w and I (intrinsic variables) are each
% column vectors with ncell components.
% Nm = number of "m"-gates, to set variance of minf.
%
% Suggested usage:
%
% cmat=[0 -1 0; 0 0 -1; -1 0 0];
% [t,y,inputs]=mlsqrnoisy3(0.2,2e4-1,cmat,1e-2,1e-2,1e3,1);

tic % record elapsed time
savefile = ['mlsqrn3_',datestr(now,30),'.mat']; % file for saving results

% a 3x3 coupling matrix, antisymmetric 
cth3 = [0 1 -1; -1 0 1; 1 -1 0]; % this leads to "pronking"

% a 3x3 coupling matrix, inhibitory and nonsymmetric
cinh3 = [0 -.5 -1; -1 0 -.5; -.5 -1 0]; % 

% a 3x3 coupling matrix, inhibitory and even less symmetric
cinh3irreg = [0 -.4 -1.1; -1 0 -.5; -.6 -.9 0]; % sum of input to each unit is -1.5

if nargin < 1, dt0 = .2; end
if nargin < 2, tmax = 1e3-1; end
if nargin < 3, cmat = cth3; end % 3x3, no coupling unless specified
global ge gi
if nargin < 4, ge = 1e-2; end % scale excitatory conductance
if nargin < 5, gi = 1e-2; end % scale inhibitory conductance
if nargin < 6, Nm = 1e3; end % number of m-gates (each cell)
%if nargin < 7, plotflag = 0; end % default = don't produce plots
if nargin < 7, plotflag = 1; end % default = do produce plots
if nargin < 8, rseed = sum(100*clock); end % seed for random number generator

inputs{1} = dt0; inputs{2}=tmax; inputs{3}=cmat; inputs{4}=ge; inputs{5}=gi; 
inputs{6} = Nm;

%% All cells have same intrinsic parameters
global epsi v0 vk vl vca
global gk gl c
global v1 v2
global v3 v4 phi gca

epsi=.001;v0=-26;vk=-84;vl=-60;vca=120;
gk=8;gl=2;c=20;
v1=-1.2;v2=18;
v3=12;v4=17.4;phi=.23;gca=4;

%% Coupling parameters (same as for ca and k -- a guess)
global ve vi
%global ge gi    % appear as inputs to function mlsqr2 -- declared above
global vsyn v6
ve = 120; vi = -84;
vsyn = -15; % threshold for synaptic input
v6 = 5;     % spread for synaptic activation

%% Number of coupled cells & type of coupling
ncell = length(cmat);
exc = sparse(sign(cmat)==1);    % boolean matrix, 1 where j->i excitatory else 0
inh = sparse(sign(cmat)==-1);   % boolean matrix, 1 where j->i inhibitory else 0

%% Initialize y = [v;w;I;Isyn]
init = '3async';
switch init
    case 'rand'
        vinit=-41.84*ones(ncell,1) + 4*randn(ncell,1);
        winit=0.002*ones(ncell,1) + .0002*randn(ncell,1);
        Iinit=30*ones(ncell,1) + 3*randn(ncell,1); % to match mean
    case 'same'
        vinit=-41.84*ones(ncell,1);
        winit=0.002*ones(ncell,1);
        Iinit=30*ones(ncell,1); % to match mlsqr.m
    case '3async'
        vinit=[-30.3211; -29.5578; -34.138];
            winit=[0.0119; 0.0035; 0.0050];
            Iinit=[40.3064; 34.3202; 38.4676];
    case '2anti'
        vinit=[-26.1210;  -37.9576];
        winit=[0.0122;    0.0032];
        Iinit=[40.5880;   35.0199];        
    case '1234'
        vinit=[-39.9070  -36.5675  -31.9199   11.2514]';
        winit=[0.0026    0.0037    0.0064    0.4055]';
        Iinit=[32.2853   36.4382   39.7143   37.9563]';
    case '4321'
        vinit=[11.2514  -31.9199  -36.5675  -39.9070]';
        winit=[0.4055    0.0064    0.0037    0.0026]';
        Iinit=[37.9563   39.7143   36.4382   32.2853]';
    case '1432'
        vinit=[-39.9070  11.2514  -31.9199   -36.5675]';
        winit=[0.0026    0.4055    0.0064     0.0037]';
        Iinit=[32.2853   37.9563   39.7143    36.4382]';        
end
inputs{7}=init;
Isyninit=zeros(ncell,1);

y0 = [vinit;winit;Iinit;Isyninit]; % a column vector of length ncell*4

%% Initialize random number generator

rand('twister',rseed);
inputs{8} = 'twister';
inputs{9} = rseed;
inputs{10} = 'mlsqrnoisy4.m';

%% Forward Euler
y = y0; v = vinit; w = winit; I = Iinit; t = 0; ts = [t]; 
while t < tmax
    if ~rem(ceil(1e6*t/tmax)/1000,1),disp([t,tmax]),end
    % scale dt to maintain arc length per time step for intrinsic dynamics
    dvdt = fv(v,w) + I/c;
    dwdt = fw(v,w);
    speed2 = (dvdt'*dvdt + dwdt'*dwdt)/ncell;   %L2 norm mean squared
    %disp([t,speed2]);
    dt = dt0/sqrt(1+speed2);    
    
    Is = Isyn(v,exc,inh);   % current from synaptic inputs
    %NOTE V IS UPDATED BEFORE fw(v,w) IS CALCULATED!!
    v = v + dt*(fvrand(v,w,Nm,ncell) + (I + Is)/c);
    w = w + dt*fw(v,w);
    I = I + dt*fI(v,w);
    y = [y,[v;w;I;Is]];
    t = t+dt;
    ts= [ts,t];
end

%% PLOTTING

if plotflag
    v = y(0*ncell+(1:ncell),:);
    w = y(1*ncell+(1:ncell),:);
    I = y(2*ncell+(1:ncell),:);
    Is= y(3*ncell+(1:ncell),:);
    t = ts; clear ts
    subplot(2,2,1)
    plot(t,v),xlabel('t'),ylabel('v')
    subplot(2,2,2)
    plot(t,w),xlabel('t'),ylabel('w')
    subplot(2,2,3)
    plot(t,I),xlabel('t'),ylabel('I')
    subplot(2,2,4)
    plot(I(1,:),v(1,:)),xlabel('I_1'),ylabel('v_1')
    shg        
end

%keyboard % for debugging...

%% END OF MAIN PROGRAM

%keyboard
inputs{end+1} = toc;
eval(['save ',savefile,' t y inputs'])

%% FUNCTIONS

function m = minf(v)
global v1 v2
m=.5*(1+tanh((v-v1)/v2));

function m = minfbinom(v,Nm,ncell)
global v1 v2
m=.5*(1+tanh((v-v1)/v2));
m=random('bino',Nm,m,ncell,1); % draw from random distribution with same mean
m=m/Nm; % rescale

function w = winf(v)
global v3 v4
w=.5*(1+tanh((v-v3)/v4));

function tau = tauw(v)
global v3 v4
tau=1./cosh((v-v3)/(2*v4));

function dvdt = fv(v,w)
global gca vca gk vk gl vl c
dvdt = (- gca*minf(v).*(v-vca)-gk*w.*(v-vk)-gl*(v-vl))/c;

function dvdtr = fvrand(v,w,Nm,ncell)
global gca vca gk vk gl vl c
dvdtr = (- gca*minfbinom(v,Nm,ncell).*(v-vca)-gk*w.*(v-vk)-gl*(v-vl))/c;

function dwdt = fw(v,w)
global phi
dwdt = phi*(winf(v)-w)./tauw(v);

function dIdt = fI(v,w)
global epsi v0
dIdt = epsi*(v0-v);

function s = se_inf(v) % NOTE INSTANTANEOUS SYNAPSE -- CHANGE TO SLOW SYNAPSE?
global vsyn v6
s = 1./(1+exp(-(v-vsyn)/v6));

function s = si_inf(v) % NOTE INSTANTANEOUS SYNAPSE -- CHANGE TO SLOW SYNAPSE?
global vsyn v6
s = 1./(1+exp(-(v-vsyn)/v6));

function I=Isyn(v,exc,inh)
global ge gi ve vi
I = -(ge*exc.*((v-ve)*se_inf(v)')...
    +gi*inh.*((v-vi)*si_inf(v)'))...
    *ones(size(v));