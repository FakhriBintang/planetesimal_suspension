% planetesimal sill rainfall: user control script
% no sticky air/space; no self gravity
% equal grid spacing
clear; %close all

RUN.ID          =  'thickcooling 1300';            % run identifier
RUN.plot        =  1;                    % switch on to plot live output
RUN.save        =  0;                    % switch on to save output files
RUN.nop         =  200;                   % output every 'nop' grid steps of transport
RUN.bnchm       =  0;                    % manufactured solution benchmark on fluid mechanics solver
RUN.diseq       =  1;                    % switch to disequilibrium approach to thermochemical evolution
%temporary
RUN.rad         =  0; %radiogenic heating

%% set model timing
NUM.yr          =  3600*24*365.25;       % seconds per year
NUM.maxstep     =  1e4;                  % maximum number of time steps
NUM.tend        =  1e8*NUM.yr;           % model stopping time [s]

% [do not modify]
NUM.dt          =  0.1*NUM.yr;          % (initial) time step [s]


%% set model domain
NUM.D           =  1000;                 % domain depth
NUM.L           =  1000;                 % domain length
NUM.N           =  150;                  % number of real x block nodes

% [do not modify]
NUM.h           =  NUM.D/NUM.N;          % spacing of x coordinates


%% set thermochemical parameters

% set initial system and component fractions
CHM.xFe0        =  0.2;                 % Fe-FeS system fraction
CHM.cFe0        =  0.22;                 % Fe-FeS fertile component fraction ([wt% S], maximum 0.35 for pure FeS
CHM.cSi0        =  0.46;                 % Si system fertile component fraction [wt% SiO2]

% set parameters
dxFe            = -0e-3;                 % amplitude of initial random perturbation to iron system
dcFe            =  0e-3;                 % amplitude of initial random perturbation to iron component
dcSi            =  0e-3;                 % amplitude of initial random perturbation to silicate component
smth            =  ((NUM.N+2)/20)^2;     % regularisation of initial random perturbation

% set phase diagram parameters
%   Fertile   ||  Refractory
CHM.TFe1    = 1000; CHM.TFe2    = 1540;  % iron system melting limits
CHM.TSi1    = 891;  CHM.TSi2    = 1839;  % silicate system melting limits
CHM.cphsSi1 = 0.4080; CHM.cphsSi2 = 0.8282;  % silicate system limits
CHM.cphsFe1 = 0     ; CHM.cphsFe2 = 0.35;  % iron system limits
CHM.perClSi = 0.5776;                      % silicate peritectic liquidus composition [wt SiO2]
CHM.perCsSi = 0.5006;                      % silicate peritectic solidus  composition [wt SiO2]
CHM.perTSi  = 1193;                      % silicate peritectic temperature
CHM.PhDgSi  = [25.9784 10.7034 0.7618 1.3523];                       % silicate phase diagram curvature factor (> 1)
CHM.perClFe = CHM.cphsFe2;               % iron peritectic liquidus composition [wt SiO2]
CHM.perCsFe = CHM.cphsFe2;               % iron peritectic solidus  composition [wt SiO2]
CHM.perTFe  = CHM.TFe1;                  % iron peritectic temperature
CHM.PhDgFe  = 5.5;                       % iron hase diagram curvature factor (> 1)
CHM.clap    = 1e-7;                      % Clapeyron slope for P-dependence of melting T [degC/Pa]

CHM.tau_r   = 1e-3*NUM.yr;                % reaction time scale [s]

% set temperature initial condition
SOL.T0      =  1600;                      % reference/top potential temperature [C]
SOL.T1      =  1600;                     % bottom potential temperature (if different from top) [C]
SOL.rT      =  NUM.D/6;                  % radius of hot plume [m]
SOL.zT      =  NUM.D*0.5;                % z-position of hot plume [m]
SOL.xT      =  NUM.L/2;                  % x-position of hot plume [m]

SOL.Ttype   = 'constant';                % set initial temperature field type

%% set material parameters
% buoyancy parameters
PHY.rhoSis      =  3300;                 % reference density solid refractory silicate [kg/m3]
PHY.rhoSil      =  2900;                 % reference density liquid refractory silicate [kg/m3]
PHY.rhoFes      =  8000;                 % reference desnity solid refractory iron [kg/m3]
PHY.rhoFel      =  7600;                 % reference desnity liquid refractory iron [kg/m3]
PHY.gCSi        =  0.50;                 % compositional expansivity silicate
PHY.gCFe        =  0.65;                 % compositional expansivity iron
PHY.aT        =  3e-5;                 % thermal expansivity silicate [1/K]
PHY.dx          =  1e-3;                 % solid grain size [m]
PHY.df          =  1e-3;                 % metal droplet size [m]
PHY.dm          =  1e-3;                 % melt film size [m]
PHY.gz          =  0.1;                  % z-gravity
PHY.gx          =  0;               	 % x-gravity
P0              = 0;

% rheology parameters
PHY.EtaSil0     =  1e2;                  % reference silicate melt viscosity [Pas]
PHY.EtaFel0     =  100;                  % reference metal melt viscosity [Pas]
PHY.EtaSol0     =  1e15;                 % reference silicate/iron crystal viscosity
Em              =  150e3;                % activation energy melt viscosity [J/mol]

AAP             =  [ 0.25, 0.25, 0.25; ...
                     0.25, 0.25, 0.25; ...
                     0.25, 0.25, 0.25; ];  % permission slopes

BBP             =  [ 0.44, 0.18, 0.38; ...
                     0.60, 0.03, 0.37; ...
                     0.70, 0.24, 0.06; ];  % permission step locations

CCP             =  [ 0.30, 0.30, 0.30; ...
                     0.60, 0.60, 0.12; ...
                     0.60, 0.12, 0.60; ];  % permission step widths

% thermochemical parameters
PHY.kTSi        =  4;                   % Thermal conductivity silicate [W/m/K]
PHY.kTFe        =  70;                  % Thermal conductivity [W/m/K]
PHY.kC          =  1e-7;                % chemical diffusivity [m^2/s]

PHY.Cp          = 1000;                 % mixture heat capacity

CHM.dEntrSi = 300;                      % silicate entropy of melting
CHM.dEntrFe = 200;                      % iron entropy of melting

PHY.Hr0         =  1e-3;                % Radiogenic heat productivity [W/m3]


%% set boundary conditions
% Temperature boundary conditions
SOL.BCTTop      = 'insulating';          % 'isothermal' or 'insulating' bottom boundaries
SOL.BCTBot      = 'insulating';          % 'isothermal' or 'insulating' bottom boundaries
SOL.BCTSides    = 'insulating';          % 'isothermal' or 'insulating' bottom boundaries

% Velocity boundary conditions: free slip = -1; no slip = 1
SOL.BCsides     = -1;                     % side boundaries
SOL.BCtop       = -1;                     % top boundary
SOL.BCbot       = -1;                     % bottom boundary


%% set solver options
% advection scheme
NUM.ADVN        = 'fromm';  % advection scheme ('fromm','first upwind','second upwind','third upwind','flxdiv')
TINY            = 1e-16;    % tiny number to safeguard [0,1] limits
NUM.CFL         = 0.2;   	% Courant number to limit physical time step
NUM.theta     	= 0.5;      % 0 = backwards Euler, 0.5 = Crank-Nicholson, 1 = Forward Euler
NUM.reltol    	= 1e-3;     % relative residual tolerance for nonlinear iterations
NUM.abstol      = 1e-6;     % absolute residual tolerance for nonlinear iterations
NUM.maxit       = 20;       % maximum iteration count
dtmax           = 0.5*NUM.yr; % maximum time step
etamin          = 1e2;      % minimum viscosity for stabilisation
etamax          = 1e15;     % maximum viscosity for stabilisation
alpha           = 0.40;     % iterative lagging parameters
nvsmooth        = 20;       % smoothing interations for the vseg boundaries. 10 rcommended for thermal boundaries, 20 for isothermal
% temporary - option for 3rd order runge-kutta
RK3 = 0;
Siseg = 0; % multiplier to liquid silicate segregation velocity, zero or one;

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

