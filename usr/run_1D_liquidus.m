% planetesimal sill rainfall: user control script
% no sticky air/space; no self gravity
% equal grid spacing
clear ; close 

RunID           =  'test';               % run identifier
plot_op         =  1;                    % switch on to plot live output
save_op         =  0;                    % switch on to save output files
nop             =  1;                   % output every 'nop' grid steps of transport
bnchm           =  0;                    % manufactured solution benchmark on fluid mechanics solver
%temporary
radheat         =  0;                    % radiogenic heating


%% set model timing
yr              =  3600*24*365.25;       % seconds per year
maxstep         =  2e4;                  % maximum number of time steps
tend            =  1e8*yr;           % model stopping time [s]

% [do not modify]
dt              =  1e-2*yr;          % (initial) time step [s]


%% set model domain
D               =  500;                  % domain depth
N               =  100;                  % number of real x/z block nodes

% [do not modify]
h               =  D/N;          % spacing of x/z  coordinates
L               =  h;


%% set thermochemical parameters

% set initial system and component fractions
xFe0            =  0.20;                 % Fe-FeS system fraction
cFe0            =  0.15;                 % Fe-FeS fertile component fraction ([wt% S], maximum 0.35 for pure FeS
cSi0            =  0.49;                 % Si system fertile component fraction [wt% SiO2]

% set parameters
dxFe            = -0e-3;                 % amplitude of initial random perturbation to iron system
dcFe            =  0e-3;                 % amplitude of initial random perturbation to iron component
dcSi            =  0e-3;                 % amplitude of initial random perturbation to silicate component
smth            =  ((N+2)/20)^2;     % regularisation of initial random perturbation

% set phase diagram parameters
%   Fertile        ||       Refractory
TFe1    = 1000;     TFe2    = 1540;   % iron system melting limits
TSi1    = 891;      TSi2    = 1839;   % silicate system melting limits
cphsSi1 = 0.4080;   cphsSi2 = 0.7882; % silicate system limits
cphsFe1 = 0     ;   cphsFe2 = 0.35;   % iron system limits
perClSi = 0.5276;                    % silicate peritectic liquidus composition [wt SiO2]
perCsSi = 0.4806;                    % silicate peritectic solidus  composition [wt SiO2]
perTSi  = 1193;                      % silicate peritectic temperature
PhDgSi  = [8.0,4.0,1.2,1.2];         % silicate phase diagram curvature factor (> 1)
perClFe = cphsFe2;               % iron peritectic liquidus composition [wt SiO2]
perCsFe = cphsFe2;               % iron peritectic solidus  composition [wt SiO2]
perTFe  = TFe1;                  % iron peritectic temperature
PhDgFe  = [8.0,4.0,1.2,1.2];         % iron hase diagram curvature factor (> 1)
clap    = 1e-7;                      % Clapeyron slope for P-dependence of melting T [degC/Pa]

% set temperature initial condition
T0      =  1600;                     % reference/top potential temperature [C]
T1      =  1600;                     % bottom potential temperature (if different from top) [C]
rT      =  D/6;                  % radius of hot plume [m]
zT      =  D*0.5;                % z-position of hot plume [m]
xT      =  L/2;                  % x-position of hot plume [m]

Ttype   = 'constant';                % set initial temperature field type


%% set material parameters
% buoyancy parameters
rhosSi      =  3300;                 % reference density solid refractory silicate [kg/m3]
rholSi      =  2900;                 % reference density liquid refractory silicate [kg/m3]
rhosFe      =  8000;                 % reference desnity solid refractory iron [kg/m3]
rholFe      =  7600;                 % reference desnity liquid refractory iron [kg/m3]
gCSi        =  0.50;                 % compositional expansivity silicate
gCFe        =  0.65;                 % compositional expansivity iron
aT          =  3e-5;                 % thermal expansivity silicate [1/K]
dx          =  1e-4;                 % solid grain size [m]
df          =  1e-4;                 % metal droplet size [m]
dm          =  1e-4;                 % melt film size [m]
gz0          =  0.1;                  % z-gravity
gx0          =  0;               	 % x-gravity

% Reference pressure
P0 = 0;

% rheology parameters
EtalSi0     =  1e2;                  % reference silicate melt viscosity [Pas]
EtalFe0     =  1e1;                  % reference metal melt viscosity [Pas]
EtaSol0     =  1e15;                 % reference silicate/iron crystal viscosity
Em          =  150e3;                % activation energy melt viscosity [J/mol]

AAP         =  [ 0.25, 0.25, 0.25; ...
                 0.25, 0.25, 0.25; ...
                 0.25, 0.25, 0.25; ];  % permission slopes

BBP         =  [ 0.44, 0.18, 0.38; ...
                 0.61, 0.01, 0.38; ...
                 0.70, 0.24, 0.06; ];  % permission step locations

CCP         =  [ 0.30, 0.30, 0.30; ...
                 0.60, 0.60, 0.12; ...
                 0.60, 0.12, 0.60; ];  % permission step widths

% thermochemical parameters
kTSi        =  3;                   % Thermal conductivity silicate [W/m/K]
kTFe        =  6;                   % Thermal conductivity [W/m/K]
kC          =  1e-7;                % chemical diffusivity [m^2/s]

Cp          = 1000;                 % mixture heat capacity

dEntrSi     = -200;                     % silicate entropy of crystallisation
dEntrFe     = -200;                     % iron-sulfide entropy of crystallisation

Hr0         =  0e-4;                % Radiogenic heat productivity [W/m3]


%% set boundary conditions
% Temperature boundary conditions
BCTTop      = 'isothermal';         % 'isothermal' or 'insulating' bottom boundaries
BCTBot      = 'insulating';         % 'isothermal' or 'insulating' bottom boundaries
BCTSides    = 'insulating';         % 'isothermal' or 'insulating' bottom boundaries

% Velocity boundary conditions: free slip = -1; no slip = 1
BCsides     = -1;                   % side boundaries
BCtop       = -1;                   % top boundary
BCbot       = -1;                   % bottom boundary


%% set solver options
% advection scheme
ADVN        =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
BCA         =  {'',''};             % boundary condition on advection (top/bot, sides)
TINY        = 1e-16;                % tiny number to safeguard [0,1] limits
theta       = 0.5;   	            % time stepping mode
CFL         = 0.25;   	            % Courant number to limit physical time step
reltol    	= 1e-6;                 % relative residual tolerance for nonlinear iterations
abstol      = 1e-9;                 % absolute residual tolerance for nonlinear iterations
maxit       = 50;                   % maximum iteration count
dtmax       = 5e-3*yr;          % maximum time step
etareg      = 1e0;                  % regularisation factor for viscosity
alpha       = 0.50;                 % iterative lagging parameters


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

run('main');

