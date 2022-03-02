% planetesimal sill rainfall: user control script
% no sticky air/space; no self gravity
% equal grid spacing
clear; close all

NN = [40,80,160];  % test increasing time steps

for nn = NN
    
RUN.ID          =  'bnchm_VP';           % run identifier
RUN.bnchm       =  1;                    % manufactured solution benchmark on fluid mechanics solver
RUN.diseq       =  1;                    % switch to disequilibrium approach to thermochemical evolution
NUM.N           =  nn;                   % number of real x block nodes


RUN.plot        =  0;                    % switch on to plot live output
RUN.save        =  0;                    % switch on to save output files
RUN.nop         =  1;                    % output every 'nop' grid steps of transport


%% set model timing
NUM.yr          =  3600*24*365.25;       % seconds per year
NUM.maxstep     =  1e4;                  % maximum number of time steps
NUM.tend        =  1e8*NUM.yr;           % model stopping time [s]

% [do not modify]
NUM.dt          =  0.1*NUM.yr;           % (initial) time step [s]


%% set model domain
NUM.D           =  1000;                 % domain depth
NUM.L           =  1000;                 % domain length

% [do not modify]
NUM.h           =  NUM.D/NUM.N;          % spacing of x coordinates


%% set thermochemical parameters

% set initial system and component fractions
CHM.xFe0        =  0.0;                  % Fe-FeS system fraction
CHM.cFe0        =  0.0;                  % Fe-FeS fertile component fraction ([wt% S], maximum 0.35 for pure FeS
CHM.cSi0        =  0.5;                  % Si system fertile component fraction [wt% SiO2]

% set parameters
dxFe            =  0;                    % amplitude of initial random perturbation to iron system
dcFe            =  0;                    % amplitude of initial random perturbation to iron component
dcSi            =  1e-3;                 % amplitude of initial random perturbation to silicate component
smth            =  ((NUM.N+2)/20)^2;     % regularisation of initial random perturbation

% set phase diagram parameters
%   Fertile   ||  Refractory
CHM.TFe1    = 1000; CHM.TFe2    = 1500;  % iron system melting limits
CHM.TSi1    = 750;  CHM.TSi2    = 1750;  % silicate system melting limits
CHM.cphsSi1 = 0.36; CHM.cphsSi2 = 0.72;  % silicate system limits
CHM.cphsFe1 = 0   ; CHM.cphsFe2 = 0.35;  % iron system limits
CHM.perClSi = 0.51;                      % silicate peritectic liquidus composition [wt SiO2]
CHM.perCsSi = 0.48;                      % silicate peritectic solidus  composition [wt SiO2]
CHM.perTSi  = 1100;                      % silicate peritectic temperature
CHM.PhDgSi  = 5.0;                       % silicate phase diagram curvature factor (> 1)
CHM.perClFe = CHM.cphsFe2;               % iron peritectic liquidus composition [wt SiO2]
CHM.perCsFe = CHM.cphsFe2;               % iron peritectic solidus  composition [wt SiO2]
CHM.perTFe  = CHM.TFe1;                  % iron peritectic temperature
CHM.PhDgFe  = 5.0;                       % iron hase diagram curvature factor (> 1)
CHM.clap    = 1e-7;                      % Clapeyron slope for P-dependence of melting T [degC/Pa]
CHM.dEntrSi = 300;                       % silicate entropy of melting
CHM.dEntrFe = 300;                       % iron entropy of melting
CHM.tau_r   = 10*DT(1);                  % reaction time scale [s]

% set temperature initial condition
SOL.T0      =  1100;                     % reference/top potential temperature [C]
SOL.T1      =  1100;                     % bottom potential temperature (if different from top) [C]
SOL.rT      =  NUM.D/6;                  % radius of hot plume [m]
SOL.zT      =  NUM.D*0.5;                % z-position of hot plume [m]
SOL.xT      =  NUM.L/2;                  % x-position of hot plume [m]

SOL.Ttype   = 'constant';                % set initial temperature field type


%% set material parameters
% buoyancy parameters
PHY.rhoSis      =  3300;                 % reference density solid refractory silicate [kg/m3]
PHY.rhoSil      =  2800;                 % reference density liquid refractory silicate [kg/m3]
PHY.rhoFes      =  8000;                 % reference desnity solid refractory iron [kg/m3]
PHY.rhoFel      =  7500;                 % reference desnity liquid refractory iron [kg/m3]
PHY.gCSi        =  0.52;                 % compositional expansivity silicate
PHY.gCFe        =  0.71;                 % compositional expansivity iron
PHY.aTSi        =  3e-5;                 % thermal expansivity silicate [1/K]
PHY.aTFe        =  1e-5;                 % thermal expansivity iron [1/K]
PHY.dx          =  1e-12;                % solid grain size [m]
PHY.df          =  1e-12;                % metal droplet size [m]
PHY.dm          =  1e-12;                % melt film size [m]
PHY.gz          =  10;                   % z-gravity
PHY.gx          =  0;               	 % x-gravity

% rheology parameters
PHY.EtaSil0     =  1e3;                  % reference silicate melt viscosity [Pas]
PHY.EtaFel0     =  1;                    % reference metal melt viscosity [Pas]
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
PHY.kC          =  1e-8;                % chemical diffusivity [m^2/s]

PHY.CpFes       =  1000;                % Volumetric heat capacity [J/kg/K]
PHY.CpFel       =  1000;                % melt heat capacity [J/kg/K]
PHY.CpSis       =  1000;                % xtal heat capacity [J/kg/K]
PHY.CpSil       =  1000;                % mvp  heat capacity [J/kg/K]

PHY.Hr0         =  0e-7;                % Radiogenic heat productivity [W/m3]


%% set boundary conditions
% Temperature boundary conditions
SOL.BCTTop      = 'insulating';          % 'isothermal' or 'insulating' bottom boundaries
SOL.BCTBot      = 'insulating';          % 'isothermal' or 'insulating' bottom boundaries
SOL.BCTSides    = 'insulating';          % 'isothermal' or 'insulating' bottom boundaries

% Velocity boundary conditions: free slip = -1; no slip = 1
SOL.BCsides     = -1;                    % side boundaries
SOL.BCtop       = -1;                    % top boundary
SOL.BCbot       = -1;                    % bottom boundary


%% set solver options
% advection scheme
NUM.ADVN        = 'fromm';  % advection scheme ('fromm','first upwind','second upwind','third upwind','flxdiv')
TINY            = 1e-16;    % tiny number to safeguard [0,1] limits
NUM.CFL         = 0.25;   	% Courant number to limit physical time step
NUM.theta     	= 0.5;      % 0 = backwards Euler, 0.5 = Crank-Nicholson, 1 = Forward Euler
NUM.reltol    	= 1e-3;     % relative residual tolerance for nonlinear iterations
NUM.abstol      = 1e-6;     % absolute residual tolerance for nonlinear iterations
NUM.maxit       = 20;       % maximum iteration count
dtmax           = 1*NUM.yr; % maximum time step
etamin          = 1e1;      % minimum viscosity for stabilisation
etamax          = 1e7;      % maximum viscosity for stabilisation
alpha           = 0.9;      % iterative lagging parameters


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


figure(17); clf;
subplot(2,3,1); imagesc(x_mms,zw_mms,- SOL.W(:,2:end-1)*hr); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. sol. $W$ [m/hr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
subplot(2,3,2); imagesc(xu_mms,z_mms,  SOL.U(2:end-1,:)*hr); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. sol. $U$ [m/hr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
subplot(2,3,3); imagesc(x_mms ,z_mms,  SOL.P(2:end-1,2:end-1)/1e3); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. sol. $P$ [kPa]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
subplot(2,3,4); imagesc(x_mms,zw_mms,-(SOL.W(:,2:end-1)-W_mms(:,2:end-1))*hr); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. err. $W$ [m/hr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
subplot(2,3,5); imagesc(xu_mms,z_mms, (SOL.U(2:end-1,:)-U_mms(2:end-1,:))*hr); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. err. $U$ [m/hr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
subplot(2,3,6); imagesc(x_mms ,z_mms, (SOL.P(2:end-1,2:end-1)-P_mms(2:end-1,2:end-1))/1e3); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. err. $P$ [kPa]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
drawnow;

% get solution error
errW = norm(SOL.W(:,2:end-1)-W_mms(:,2:end-1),2)./norm(W_mms(:,2:end-1),2);
errU = norm(SOL.U(2:end-1,:)-U_mms(2:end-1,:),2)./norm(U_mms(2:end-1,:),2);
errP = norm(SOL.P(2:end-1,2:end-1)-P_mms(2:end-1,2:end-1),2)./norm(P_mms(2:end-1,2:end-1),2);

% plot error convergence
figure(18); 
p1 = loglog(h,errW,'r+','MarkerSize',8,'LineWidth',2); axis xy tight; hold on; box on; 
p2 = loglog(h,errU,'g+','MarkerSize',8,'LineWidth',2); axis xy tight; hold on; box on; 
p3 = loglog(h,errP,'b+','MarkerSize',8,'LineWidth',2); axis xy tight; hold on; box on; 
set(gca,'TicklabelInterpreter','latex','FontSize',12)
xlabel('grid spacing [m]','Interpreter','latex')
ylabel('numerical error [1]','Interpreter','latex')
set(gca,'TicklabelInterpreter','latex')
title('Numerical convergence in space','Interpreter','latex','FontSize',20)
    
if nn == NN(1)
    p4 = loglog(L./NN,mean([errW,errU,errP]).*(NN(1)./NN).^2,'k-','LineWidth',2);  % plot linear trend for comparison
end
if nn == NN(end)
    legend([p1,p2,p3,p4],{'error W','error U','error P','quadratic'},'Interpreter','latex','box','on','location','southeast')
end

% plot error convergence
figure(19);
DOFS = (NN+2).*(NN+2) + 2.*(NN+1).*(NN+2);
dofs = (nn+2).*(nn+2) + 2.*(nn+1).*(nn+2);
p5 = loglog(dofs,solvetime,'r+','MarkerSize',8,'LineWidth',2); axis xy tight; hold on; box on; 
set(gca,'TicklabelInterpreter','latex','FontSize',12)
xlabel('\# dofs [1]','Interpreter','latex','FontSize',16)
ylabel('time to solution [s]','Interpreter','latex','FontSize',16)
title('Scaling of direct solver','Interpreter','latex','FontSize',20)

if nn == NN(1)
    p6 = loglog(DOFS,0.95*solvetime*(DOFS./DOFS(1)).^1,'k-','LineWidth',2);  % plot linear trend for comparison
end
if nn == NN(end)
    legend([p5,p6],{'time to solution','linear'},'Interpreter','latex','box','on','location','southeast')
end

end

name = [outdir,'/',RUN.ID,'/',RUN.ID,'_bnchm'];
print(fh18,name,'-dpng','-r300','-opengl');

name = [outdir,'/',RUN.ID,'/',RUN.ID,'_sclng'];
print(fh19,name,'-dpng','-r300','-opengl');


