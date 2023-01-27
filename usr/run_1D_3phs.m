% planetesimal sill rainfall: user control script
% no sticky air/space; no self gravity
% equal grid spacing
clear all; close all

RUN.ID          =  'test';               % run identifier
RUN.plot        =  1;                    % switch on to plot live output
RUN.save        =  0;                    % switch on to save output files
RUN.nop         =  10;                   % output every 'nop' grid steps of transport
RUN.bnchm       =  0;                    % manufactured solution benchmark on fluid mechanics solver
RUN.diseq       =  1;                    % switch to disequilibrium approach to thermochemical evolution
%temporary
RUN.rad         =  0;                    % radiogenic heating


%% set model timing
NUM.yr          =  3600*24*365.25;       % seconds per year
NUM.maxstep     =  2e4;                  % maximum number of time steps
NUM.tend        =  1e8*NUM.yr;           % model stopping time [s]

% [do not modify]
NUM.dt          =  1e-2*NUM.yr;          % (initial) time step [s]


%% set model domain
NUM.D           =  500;                  % domain depth
NUM.N           =  100;                  % number of real x/z block nodes

% [do not modify]
NUM.h           =  NUM.D/NUM.N;          % spacing of x/z  coordinates
NUM.L           =  NUM.h;


%% set thermochemical parameters

% set initial system and component fractions
CHM.xFe0        =  0.20;                 % Fe-FeS system fraction
CHM.cFe0        =  0.15;                 % Fe-FeS fertile component fraction ([wt% S], maximum 0.35 for pure FeS
CHM.cSi0        =  0.46;                 % Si system fertile component fraction [wt% SiO2]

% set parameters
dxFe            = -0e-3;                 % amplitude of initial random perturbation to iron system
dcFe            =  0e-3;                 % amplitude of initial random perturbation to iron component
dcSi            =  0e-3;                 % amplitude of initial random perturbation to silicate component
smth            =  ((NUM.N+2)/20)^2;     % regularisation of initial random perturbation

% set phase diagram parameters
%   Fertile        ||       Refractory
CHM.TFe1    = 1000;     CHM.TFe2    = 1540;   % iron system melting limits
CHM.TSi1    = 891;      CHM.TSi2    = 1839;   % silicate system melting limits
CHM.cphsSi1 = 0.4080;   CHM.cphsSi2 = 0.7882; % silicate system limits
CHM.cphsFe1 = 0     ;   CHM.cphsFe2 = 0.35;   % iron system limits
CHM.perClSi = 0.5276;                    % silicate peritectic liquidus composition [wt SiO2]
CHM.perCsSi = 0.4806;                    % silicate peritectic solidus  composition [wt SiO2]
CHM.perTSi  = 1193;                      % silicate peritectic temperature
CHM.PhDgSi  = [8.0,4.0,1.2,1.2];         % silicate phase diagram curvature factor (> 1)
CHM.perClFe = CHM.cphsFe2;               % iron peritectic liquidus composition [wt SiO2]
CHM.perCsFe = CHM.cphsFe2;               % iron peritectic solidus  composition [wt SiO2]
CHM.perTFe  = CHM.TFe1;                  % iron peritectic temperature
CHM.PhDgFe  = [8.0,4.0,1.2,1.2];         % iron hase diagram curvature factor (> 1)
CHM.clap    = 1e-7;                      % Clapeyron slope for P-dependence of melting T [degC/Pa]

% set temperature initial condition
SOL.T0      =  1600;                     % reference/top potential temperature [C]
SOL.T1      =  1000;                     % bottom potential temperature (if different from top) [C]
SOL.rT      =  NUM.D/6;                  % radius of hot plume [m]
SOL.zT      =  NUM.D*0.5;                % z-position of hot plume [m]
SOL.xT      =  NUM.L/2;                  % x-position of hot plume [m]

SOL.Ttype   = 'constant';                % set initial temperature field type


%% set material parameters
% buoyancy parameters
PHY.rhosSi      =  3300;                 % reference density solid refractory silicate [kg/m3]
PHY.rholSi      =  2900;                 % reference density liquid refractory silicate [kg/m3]
PHY.rhosFe      =  8000;                 % reference desnity solid refractory iron [kg/m3]
PHY.rholFe      =  7600;                 % reference desnity liquid refractory iron [kg/m3]
PHY.gCSi        =  0.50;                 % compositional expansivity silicate
PHY.gCFe        =  0.65;                 % compositional expansivity iron
PHY.aT          =  3e-5;                 % thermal expansivity silicate [1/K]
PHY.dx          =  1e-4;                 % solid grain size [m]
PHY.df          =  1e-4;                 % metal droplet size [m]
PHY.dm          =  1e-4;                 % melt film size [m]
PHY.gz          =  0.1;                  % z-gravity
PHY.gx          =  0;               	 % x-gravity

% Reference pressure
P0 = 0;

% rheology parameters
PHY.EtalSi0     =  1e2;                  % reference silicate melt viscosity [Pas]
PHY.EtalFe0     =  1e1;                  % reference metal melt viscosity [Pas]
PHY.EtaSol0     =  1e15;                 % reference silicate/iron crystal viscosity
Em              =  150e3;                % activation energy melt viscosity [J/mol]

AAP             =  [ 0.25, 0.25, 0.25; ...
                     0.25, 0.25, 0.25; ...
                     0.25, 0.25, 0.25; ];  % permission slopes

BBP             =  [ 0.44, 0.18, 0.38; ...
                     0.61, 0.01, 0.38; ...
                     0.70, 0.24, 0.06; ];  % permission step locations

CCP             =  [ 0.30, 0.30, 0.30; ...
                     0.60, 0.60, 0.12; ...
                     0.60, 0.12, 0.60; ];  % permission step widths

% thermochemical parameters
PHY.kTSi        =  3;                   % Thermal conductivity silicate [W/m/K]
PHY.kTFe        =  6;                   % Thermal conductivity [W/m/K]
PHY.kC          =  1e-7;                % chemical diffusivity [m^2/s]

PHY.Cp          = 1000;                 % mixture heat capacity

CHM.dEntrSi = -200;                     % silicate entropy of crystallisation
CHM.dEntrFe = -200;                     % iron-sulfide entropy of crystallisation

PHY.Hr0         =  0e-4;                % Radiogenic heat productivity [W/m3]


%% set boundary conditions
% Temperature boundary conditions
SOL.BCTTop      = 'isothermal';         % 'isothermal' or 'insulating' bottom boundaries
SOL.BCTBot      = 'insulating';         % 'isothermal' or 'insulating' bottom boundaries
SOL.BCTSides    = 'insulating';         % 'isothermal' or 'insulating' bottom boundaries

% Velocity boundary conditions: free slip = -1; no slip = 1
SOL.BCsides     = -1;                   % side boundaries
SOL.BCtop       = -1;                   % top boundary
SOL.BCbot       = -1;                   % bottom boundary


%% set solver options
% advection scheme
ADVN            =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
BCA             =  {'',''};             % boundary condition on advection (top/bot, sides)
TINY            = 1e-16;                % tiny number to safeguard [0,1] limits
NUM.theta       = 0.5;   	            % time stepping mode
NUM.CFL         = 0.25;   	            % Courant number to limit physical time step
NUM.reltol    	= 1e-4;                 % relative residual tolerance for nonlinear iterations
NUM.abstol      = 1e-6;                 % absolute residual tolerance for nonlinear iterations
NUM.maxit       = 50;                   % maximum iteration count
dtmax           = 5e-3*NUM.yr;          % maximum time step
etareg          = 1e0;                  % regularisation factor for viscosity
alpha           = 0.50;                 % iterative lagging parameters


%% start model
% create output directory
[~,systemname] = system('hostname');
systemname(end) = [];
switch systemname
    case 'Horatio'
        outpath = ['/media/43TB_RAID_Array/fbintang/test_out/out/', RUN.ID];
        if ~exist(outpath, 'dir'); mkdir(outpath); end
    otherwise
        outpath = ['../out/',RUN.ID];
        if ~exist(outpath, 'dir'); mkdir(outpath); end
end

% add path to source directory
addpath('../src')
addpath('../src/cbrewer/')

% use color brewer to create colormaps
cm1 =        cbrewer('seq','YlOrRd',30) ; % sequential colour map
cm2 = flipud(cbrewer('div','RdBu'  ,30)); % divergent colour map

run('main');

