% planetesimal sill rainfall: user control script
% no sticky air/space; no self gravity
% equal grid spacing
clear ; close all

RunID           =  '2D_heating_5';               % run identifier
plot_op         =  1;                    % switch on to plot live output
save_op         =  1;                    % switch on to save output files
nop             =  1000;                   % output every 'nop' grid steps of transport
bnchm           =  0;                    % manufactured solution benchmark on fluid mechanics solver

%% set model timing
yr              =  3600*24*365.25;       % seconds per year
maxstep         =  100000;                  % maximum number of time steps
tend            =  1e10*yr;           % model stopping time [s]

% [do not modify]
dt              =  1e-3*yr;          % (initial) time step [s]

%% set model domain
D               =  100000;                  % domain depth
Nz              =  150;                  % number of real x/z block nodes
Nx              = 20;
% [do not modify]
h               =  D/Nz;          % spacing of x/z  coordinates
L               =  h*Nx;

%% set thermochemical parameters

% set initial system and component fractions
xFe0            =  0.2;                     % Fe-FeS system fraction
cFe0            =  0.15;                    % Fe-FeS fertile component fraction ([wt% S], maximum 0.35 for pure FeS
cSi0            =  0.47;                    % Si system fertile component fraction [wt% SiO2]

% set parameters
dxFe            = -0.0e-3;                  % amplitude of initial random perturbation to iron system
dcFe            =  0e-3;                    % amplitude of initial random perturbation to iron component
dcSi            =  0e-3;                    % amplitude of initial random perturbation to silicate component
smth            =  ((Nz+2)/20)^2;            % regularisation of initial random perturbation

% set phase diagram parameters
%      Fertile        ||       Refractory
TFe1    = 1000+273.15;     TFe2    = 1540+273.15;   % iron system melting limits [k]
TSi1    = 1193+273.15;     TSi2    = 1839+273.15;   % silicate system melting limits [k]
cSimin  = 0.4080;                                   % reference cSi
cphsSi1 = 0;        cphsSi2 = 0.5276-cSimin;        % silicate system limits
cphsFe1 = 0     ;   cphsFe2 = 0.35;                 % iron system limits
perclSi = cphsSi2;                                  % silicate peritectic liquidus composition [wt SiO2]
percsSi = cphsSi2;                                  % silicate peritectic solidus  composition [wt SiO2]
perTSi  = TSi1;                                     % silicate peritectic temperature
PhDgSi  = [8.0,4.0,1.2,1.2];                        % silicate phase diagram curvature factor (> 1)
perclFe = cphsFe2;                                  % iron peritectic liquidus composition [wt SiO2]
percsFe = cphsFe2;                                  % iron peritectic solidus  composition [wt SiO2]
perTFe  = TFe1;                                     % iron peritectic temperature
PhDgFe  = [8.0,4.0,1.2,1.2];                        % iron hase diagram curvature factor (> 1)
clap    = 1e-7;                                     % Clapeyron slope for P-dependence of melting T [degC/Pa]

% set temperature initial condition
T0      =  0+273.15;                             % reference/top potential temperature [k]
Ttop0   =  273.15;                                      % isothermal top reference temperature 
T1      =  0+273.15;                             % bottom potential temperature (if different from top) [k]
Tbot0   =  T1;                                      % isothermal bottom reference temperature 
rT      =  D/6;                                     % radius of hot plume [m]
zT      =  D*0.5;                                   % z-position of hot plume [m]
xT      =  L/2;                                     % x-position of hot plume [m]

Ttype   = 'constant';                               % set initial temperature field type

%% set material parameters
% buoyancy parameters
rhosSi0     =  3300;                 % reference density solid refractory silicate [kg/m3]
rholSi0     =  2900;                 % reference density liquid refractory silicate [kg/m3]
rhosFe0     =  8000;                 % reference desnity solid refractory iron [kg/m3]
rholFe0     =  7600;                 % reference desnity liquid refractory iron [kg/m3]
gCSi        =  0.50;                 % compositional expansivity silicate
gCFe        =  0.65;                 % compositional expansivity iron
aT          =  3e-5;                 % thermal expansivity silicate [1/K]
dx          =  1e-4;                 % solid grain size [m]
df          =  1e-4;                 % metal droplet size [m]
dm          =  1e-4;                 % melt film size [m]
gz0         =  0.1;                  % z-gravity
gx0         =  0;               	 % x-gravity

% Reference pressure
P0 = 0;

% rheology parameters
EtalSi0     =  1e2;                     % reference silicate melt viscosity [Pas]
EtalFe0     =  1e1;                     % reference metal melt viscosity [Pas]
EtaSol0     =  1e15;                    % reference silicate/iron crystal viscosity
Em          =  150e3;                   % activation energy melt viscosity [J/mol]

AAP         =  [ 0.25, 0.25, 0.25; ...
                 0.25, 0.25, 0.25; ...
                 0.25, 0.25, 0.25; ];   % permission slopes

BBP         =  [ 0.44, 0.18, 0.38; ...
                 0.61, 0.01, 0.38; ...
                 0.70, 0.24, 0.06; ];   % permission step locations

CCP         =  [ 0.30, 0.30, 0.30; ...
                 0.60, 0.60, 0.12; ...
                 0.60, 0.12, 0.60; ];   % permission step widths

% thermochemical parameters
kTSi        =  3;                       % Thermal conductivity silicate [W/m/K]
kTFe        =  6;                       % Thermal conductivity [W/m/K]
kC          =  1e-7;                    % chemical diffusivity [m^2/s]

Cp          = 1000;                     % mixture heat capacity

dEntrSi     = -200;                     % silicate entropy of crystallisation
dEntrFe     = -200;                     % iron-sulfide entropy of crystallisation

%% set heating parameters (if turned on)
radheat = 1;
Hr0         =  0e-4;                    % constant Radiogenic heat productivity [W/kg]
% Dynamic radiogenic heating rate
if radheat
    t_form      = 0.0*yr*1e6;     % planetesimal formation time after CAI, recommend no more than 2 half lives
    mr_Al       = 27;           % atomic mass of Al [g/mol]
    AV          = 6.022e23;     % avogadros number [mol^{-1}]
    nAl_C       = 2.62e23;      % chondritic abundance of Al [kg^{-1}]
    Al26_27     = 5.25e-5;      % cannonical initial ratio of 26Al/27Al at CAI formation
    t_halfAl    = 717000*yr; % half life of 26Al
    EAl         = 5e-13;        % decay energy
end

%% set boundary conditions
% Temperature boundary conditions
BCTTop      = 'isothermal';             % 'isothermal', 'insulating', or 'flux' bottom boundaries
BCTBot      = 'insulating';             % 'isothermal', 'insulating', or 'flux' bottom boundaries
BCTSides    = 'insulating';             % 'isothermal' or 'insulating' bottom boundaries

% Velocity boundary conditions: free slip = -1; no slip = 1
BCsides     = -1;                       % side boundaries
BCtop       = -1;                       % top boundary
BCbot       = -1;                       % bottom boundary
switch BCTTop
    case'flux'
    qT0 = -10;
end

%% set solver options
% advection scheme
ADVN        =  'weno5';                 % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
BCA         =  {'',''};                 % boundary condition on advection (top/bot, sides)
TINY        = 1e-16;                    % tiny number to safeguard [0,1] limits
lambda      = 0.5;   	                % iterative lagging for phase fraction
reltol    	= 1e-3;                     % relative residual tolerance for nonlinear iterations
abstol      = 1e-6;                     % absolute residual tolerance for nonlinear iterations
maxit       = 20;                       % maximum iteration count
CFL         = 0.25;                     % (physical) time stepping courant number (multiplies stable step) [0,1]
dtmax       = 5e2*yr;                   % maximum time step
etareg      = 1e5;                      % regularisation factor for viscosity
TINT        =  'bd3i';                  % time integration scheme ('bwei','cnsi','bd3i','bd3s')

%% start model
% create output directory
[~,systemname]  = system('hostname');
systemname(end) = [];
switch systemname
    case 'Horatio'
        outpath = ['/media/43TB_RAID_Array/fbintang/test_out/out/', RunID];
        if ~exist(outpath, 'dir'); mkdir(outpath); end
    otherwise
        outpath = ['../out/',RunID];
        if ~exist(outpath, 'dir'); mkdir(outpath); end
end

% add path to source directory
addpath('../src')
addpath('../src/cbrewer/')

% use color brewer to create colormaps
cm1 =        cbrewer('seq','YlOrRd',30) ; % sequential colour map
cm2 = flipud(cbrewer('div','RdBu'  ,30)); % divergent colour map
load ocean.mat;

% print run header
fprintf(1,'\n\n************************************************************\n');
fprintf(1,    '*****  planetesimal  |  %s  |  %s  *****\n'         ,RunID,datetime);
fprintf(1,    '************************************************************\n\n');

initialise;
run('main');


