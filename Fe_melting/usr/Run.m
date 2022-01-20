% planetesimal sill rainfall: user control script
% no sticky air/space; no self gravity
% equal grid spacing
clear; close all

RUN.ID          =  'fully molten';              % run identifier
RUN.plot        =  1;                   % switch on to plot live output
RUN.save        =  1;                   % switch on to save output files
RUN.nop         =  10;                  % output every 'nop' grid steps of transport
RUN.nup         =  1;                   % update every 'nup' grid steps of transport

%% set model timing
NUM.yr          =  3600*24*365.25;      % seconds per year
NUM.maxstep     =  1e4;                 % maximum number of time steps
NUM.tend        =  1e8*NUM.yr;          % model stopping time [s]

% [do not modify]
NUM.dt          =  1e3*NUM.yr;          % (initial) time step [s]

%% set model domain
NUM.D           =  100;                 % domain depth
NUM.L           =  100;                 % domain length
NUM.N           =  100;                 % number of real x block nodes

% [do not modify]
NUM.h           =  NUM.D/NUM.N;         % spacing of x coordinates

%% set thermochemical parameters
% set initial fractions
CHM.xFe0            =  0.3;                 % Fe-FeS system fraction
CHM.cFe0            =  0.2;                % Fe-FeS fertile component fraction ([Wt% S], maximum 0.35 for pure FeS
CHM.cSi0            =  0.52;                % Si system fertile component fraction [Wt% SiO2]

% set parameters
dc              =  0.001;              % amplitude of random noise

% set phase diagram parameters
%   Fertile   ||  Refractory
CHM.TFe1    = 1000; CHM.TFe2    = 1500;  % iron system melting limits
CHM.TSi1    = 750;  CHM.TSi2    = 1750;  % Silicate system melting limits
CHM.cphsSi1 = 0.36; CHM.cphsSi2 = 0.72;  % silicate system limits
CHM.cphsFe1 = 0   ; CHM.cphsFe2 = 0.35;  % iron system limits
CHM.perClSi  =  0.54;                % peritectic liquidus composition [wt SiO2]
CHM.perCsSi  =  0.50;                % peritectic solidus  composition [wt SiO2]
CHM.perTSi   =  1100;                 % peritectic temperature
CHM.PhDgSi   =  4.0;                 % Phase diagram curvature factor (> 1)
CHM.perClFe  =  0.35;                % peritectic liquidus composition [wt SiO2]
CHM.perCsFe  =  0.35;                % peritectic solidus  composition [wt SiO2]
% CHM.perClFe  =  0;                % peritectic liquidus composition [wt S]
% CHM.perCsFe  =  0;                % peritectic solidus  composition [wt S]
CHM.perTFe   =  1000;                 % peritectic temperature
CHM.PhDgFe   =  2.0;                 % Phase diagram curvature factor (> 1)
CHM.clap     =  0;
% CHM.clap     =  1e-7;                % Clapeyron slope for P-dependence of melting T [degC/Pa]
CHM.dEntrSi  = 300;                  % entropy of fusion
CHM.dEntrFe  = 300;
% CHM.dEntrSi  = 0.7;                  % entropy of fusion
% CHM.dEntrFe  = 0.5;

SOL.T0          =  250;                % reference/top potential temperature [C]
SOL.T1          =  1600;            	% bottom potential temperature (if different from top) [C]
SOL.dT          =  1;                   % temperature perturbation amplitude [C]
SOL.rT          =  NUM.D/6;             % radius of hot plume [m]
SOL.zT          =  NUM.D*0.5;           % z-position of hot plume [m]
SOL.xT          =  NUM.L/2;             % x-position of hot plume [m]

SOL.Ttype       = 'constant';           % constant ambient background temperature
Topbound        = 'surface';
% SOL.Ttype       = 'gaussian';             % linear temperaure profile between top and bottom

%% set material parameters
smth            =  ((NUM.N+2)/20)^2;    % regularisation of initial random perturbation
PHY.drho        =  250;                 % solid-liquid density contrast
PHY.rhoSis      =  3200;     
PHY.rhoSif      =  2600;     
PHY.rhoFes      =  8000;                % reference solid/pure refractory Fe density [kg/m3]
PHY.rhoFef      =  6000;                % pure fertile FeS density     
PHY.rhoSil      =  PHY.rhoSis-PHY.drho; % reference liquid silicate density   
PHY.rhoFel      =  PHY.rhoFes-PHY.drho; % reference iron density
PHY.gammaSi     =  (PHY.rhoSis - PHY.rhoSif)/(CHM.cphsSi2-CHM.cphsSi1)/PHY.rhoSis; % assume zero density contrast on silicates for now
PHY.gammaFe     =  (PHY.rhoFes - PHY.rhoFef)/(CHM.cphsFe2-CHM.cphsFe1)/PHY.rhoFes;
PHY.Eta0        =  1e4;                 % reference viscosity [Pas]
PHY.aT0         =  3e-5;                % thermal expansivity [1/K]
PHY.kTSi        =  4;                   % Thermal conductivity silicate [W/m/K]
PHY.kTFe        =  70;                   % Thermal conductivity [W/m/K]
PHY.kC          =  1e-9;                % chemical diffusivity [m^2/s]

PHY.CpFes       =  1000;                % Volumetric heat capacity [J/kg/K]
PHY.CpFel       =  1000;                % melt heat capacity [J/kg/K]
PHY.CpSis       =  1000;                % xtal heat capacity [J/kg/K]
PHY.CpSil       =  1000;                % mvp  heat capacity [J/kg/K]

PHY.Hr0         =  1e-4;                % Radiogenic heat productivity [W/m3]
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
ADVN = 'fromm';
% NUM.AdvnScheme  = 'first upwind'
% NUM.AdvnScheme  = 'second upwind';
% NUM.AdvnScheme  = 'third upwind';
% NUM.AdvnScheme  = 'flxdiv'
TINY  = 1e-16;
NUM.CFL         = 0.5;   	% Courant number to limit physical time step
NUM.theta     	= 0.5;      % 0 = backwards Euler, 0.5 = Crank-Nicholson, 1 = Forward Euler
NUM.restol    	= 1e-4;     % residual tolerance for nonlinear iterations
NUM.abstol      = 1e-4;
NUM.maxit       = 15;
alpha           = 0.5;
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
run('initialise_melting');

% check thermal rayleigh number
% PHY.RhoT0       = (PHY.RhoSi0*(1-PHY.phi0)+PHY.RhoFe0.*PHY.phi0);
% Ra              =  PHY.RhoT0 * PHY.aT0*(SOL.T1-SOL.T0)*(NUM.D^3)...
%                   *gz/PHY.Eta0/(PHY.kT0/PHY.RhoT0/PHY.Cp0);
% if Ra < 1e3
%     disp('WARNING: Rayleigh number too low, no free convection')
% end

% compare thermal and phase buoyancy terms
% dRhoT           = PHY.aT0.*(max(SOL.T(:))-min(SOL.T(:))).*PHY.RhoT0; 
% dRhoC           = (max(SOL.phi(:))-min(SOL.phi(:))).*(PHY.RhoFe0-PHY.RhoSi0);

% if dRhoT>dRhoC; disp('Thermal buoyancy driven run');
% else            disp('Phase buoyancy driven run');
% end

run('main');

