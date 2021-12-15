% planetesimal sill rainfall: user control script
% no sticky air/space; no self gravity
% equal grid spacing
clear; close all

RUN.ID          =  'demo';              % run identifier
RUN.plot        =  1;                   % switch on to plot live output
RUN.save        =  0;                   % switch on to save output files
RUN.nop         =  1;                   % output every 'nop' grid steps of transport
RUN.nup         =  1;                   % update every 'nup' grid steps of transport

%% set model timing
NUM.yr          =  3600*24*365.25;      % seconds per year
NUM.maxstep     =  1e4;                 % maximum number of time steps
NUM.tend        =  1e8*NUM.yr;          % model stopping time [s]

% [do not modify]
NUM.dt          =  1e3*NUM.yr;          % (initial) time step [s]

%% set model domain
NUM.D           =  10000;               % domain depth
NUM.L           =  10000;               % domain length
NUM.N           =  100;                 % number of real x block nodes

% [do not modify]
NUM.h           =  NUM.D/NUM.N;         % spacing of x coordinates

%% set thermochemical parameters
CHM.Si0         =  0.9;                 % initial average bulk silicate wt%
CHM.FeS0        =  0.1;                 % initial average bulk FeS wt%
CHM.xS          =  0.2;                 % wt% of S in FeS
dc              =  0.001;              % amplitude of random noise

SOL.Ta          =  273;                 % reference air/space temperature
SOL.T0          =  1200;                % reference/top potential temperature [C]
SOL.T1          =  1200;            	% bottom potential temperature (if different from top) [C]
SOL.dT          =  100;                 % temperature perturbation amplitude [C]
SOL.rT          =  NUM.D/6;             % radius of hot plume [m]
SOL.zT          =  NUM.D*0.5;           % z-position of hot plume [m]
SOL.xT          =  NUM.L/2;             % x-position of hot plume [m]

SOL.Ttype       = 'constant';           % constant ambient background temperature
SOL.Ttype       = 'linear';             % linear temperaure profile between top and bottom

%% set material parameters
smth            =  ((NUM.N+2)/25)^2;    % regularisation of initial random perturbation

PHY.RhoSi0      =  2500;     PHY.RhoFe0     = 7000;            % reference density [kg/m3]
PHY.Eta0        =  1e5;                 % reference viscosity [Pas]
PHY.aT0         =  3e-5;                % thermal expansivity [1/K]
PHY.kT0         =  4;                   % Thermal conductivity [W/m/K]
PHY.Cp0         =  1000;                % Volumetric heat capacity [J/kg/K]
PHY.Hr0         =  1e-8;                % Radiogenic heat productivity [W/m3]
PHY.gz          =  0.1;                 % z-gravity 
gz              = PHY.gz;
PHY.gx          =  0;               	% x-gravity

%% droplet parameters
PHY.x1          =  0.1;                 % component at the bottom boundary
PHY.x0          =  0.1;                 % component at the top boundary
PHY.phi0        =  0.1;                 % initial droplet fraction
PHY.d           =  0.01;                % droplet radius
SOL.philim      =  1e-4;                % limit liquid fraction for numerical stability
zlay            =  0.1;                 % layer thickness of top          
 
%% set boundary conditions
% Temperature boundary conditions
SOL.BCTTop     = 'isothermal';          % 'isothermal' or 'insulating' bottom boundaries
SOL.BCTBot     = 'isothermal';          % 'isothermal' or 'insulating' bottom boundaries
SOL.BCTSides   = 'insulating';          % 'isothermal' or 'insulating' bottom boundaries

% Velocity boundary conditions: free slip = -1; no slip = 1
SOL.BCleft      = -1;                   % left side boundary
SOL.BCright     = -1;                   % right side boundary
SOL.BCtop       = 1;                    % top boundary
SOL.BCbot       = 1;                    % bottom boundary
% segregation boundary conditions
SOL.BCSeg       = 1;                    % 1 = free-slip walls; 2 = no slip walls

%% set solver options
% advection scheme
NUM.AdvnScheme  = 'fromm';
% NUM.AdvnScheme  = 'first upwind'
% NUM.AdvnScheme  = 'second upwind';
% NUM.AdvnScheme  = 'third upwind';
% NUM.AdvnScheme  = 'flxdiv'

NUM.CFL         = 0.5;   	% Courant number to limit physical time step
NUM.theta     	= 0.5;      % 0 = backwards Euler, 0.5 = Crank-Nicholson, 1 = Forward Euler
NUM.restol    	= 1e-3;     % residual tolerance for nonlinear iterations
NUM.abstol      = 1e-6;
NUM.maxit       = 20;
NUM.cstab     	= 1e-6;     % stabilising coefficient for P-diagonal

%% start model
% create output directory
if ~exist(['../out/',RUN.ID],'dir'); mkdir(['../out/',RUN.ID]); end

% add path to source directory
addpath('../src')
addpath('../src/cbrewer/')

% use color brewer to create colormaps
cm1 =        cbrewer('seq','YlOrRd',30) ; % sequential colour map
cm2 = flipud(cbrewer('div','RdBu'  ,30)); % divergent colour map

% run model
run('initialise_nomelt');

% check thermal rayleigh number
PHY.RhoT0       = (PHY.RhoSi0*(1-PHY.phi0)+PHY.RhoFe0.*PHY.phi0);
Ra              =  PHY.RhoT0 * PHY.aT0*(SOL.T1-SOL.T0)*(NUM.D^3)...
                  *gz/PHY.Eta0/(PHY.kT0/PHY.RhoT0/PHY.Cp0);
if Ra < 1e3
    disp('WARNING: Rayleigh number too low, no free convection')
end

% compare thermal and phase buoyancy terms
dRhoT           = PHY.aT0.*(max(SOL.T(:))-min(SOL.T(:))).*PHY.RhoT0; 
dRhoC           = (max(SOL.phi(:))-min(SOL.phi(:))).*(PHY.RhoFe0-PHY.RhoSi0);

if dRhoT>dRhoC; disp('Thermal buoyancy driven run');
else            disp('Phase buoyancy driven run');
end

run('main');

