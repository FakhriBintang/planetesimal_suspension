% planetesimal sill rainfall: user control script
% no sticky air/space; no self gravity
% equal grid spacing
clear all; close all
RunID           =  '2D_bnchm_dt_bd3i_4phs';               % run identifier
% create output directory
[~,systemname]  = system('hostname');
systemname(end) = [];
switch systemname
    case 'Horatio'
        outpath = ['/media/43TB_RAID_Array/fbintang/test_out/out/benchmark_dt/',RunID];
        if ~exist(outpath, 'dir'); mkdir(outpath); end
    otherwise
        outpath = ['../out/benchmark_dt/', RunID];
        if ~exist(outpath, 'dir'); mkdir(outpath); end
end

% add path to source directory
addpath('../src')
addpath('../src/cbrewer/')

% use color brewer to create colormaps
cm1 =        cbrewer('seq','YlOrRd',30) ; % sequential colour map
cm2 = flipud(cbrewer('div','RdBu'  ,30)); % divergent colour map

plot_op         =  1;                    % switch on to plot live output
save_op         =  0;                    % switch on to save output files
nop             =  100; %D/dt;                   % output every 'nop' grid steps of transport
bnchm           =  0;                    % manufactured solution benchmark on fluid mechanics solver
%temporary
radheat         =  0;                    % radiogenic heating

%% set model domain
N               =  100;                  % number of real nodes
D               =  100;                  % domain depth
% [do not modify]
h               =  D/N;          % spacing of x/z  coordinates
L               =  D;

%% set thermochemical parameters
% set initial system and component fractions
xFe0            =  0.5;                 % Fe-FeS system fraction
cFe0            =  0.15;                 % Fe-FeS fertile component fraction ([wt% S], maximum 0.35 for pure FeS
cSi0            =  0.42;                 % Si system fertile component fraction [wt% SiO2]

% set parameters
dxFe            =  0.1;                 % amplitude of initial random perturbation to iron system
dcFe            =  0e-3;                 % amplitude of initial random perturbation to iron component
dcSi            =  0e-3;                 % amplitude of initial random perturbation to silicate component

% set phase diagram parameters
%   Fertile        ||       Refractory
TFe1    = 1000;     TFe2    = 1540;   % iron system melting limits
TSi1    = 891;      TSi2    = 1839;   % silicate system melting limits
cSimin  = 0.4080;                    % reference cSi (In testing)
cphsSi1 = 0;        cphsSi2 = 0.7882-cSimin; % silicate system limits
cphsFe1 = 0     ;   cphsFe2 = 0.35;   % iron system limits
perclSi = 0.5276-cSimin;                    % silicate peritectic liquidus composition [wt SiO2]
percsSi = 0.4806-cSimin;                    % silicate peritectic solidus  composition [wt SiO2]
perTSi  = 1193;                      % silicate peritectic temperature
PhDgSi  = [8.0,4.0,1.2,1.2];         % silicate phase diagram curvature factor (> 1)
perclFe = cphsFe2;               % iron peritectic liquidus composition [wt SiO2]
percsFe = cphsFe2;               % iron peritectic solidus  composition [wt SiO2]
perTFe  = TFe1;                  % iron peritectic temperature
PhDgFe  = [8.0,4.0,1.2,1.2];         % iron hase diagram curvature factor (> 1)
clap    = 1e-7;                      % Clapeyron slope for P-dependence of melting T [degC/Pa]

% set temperature initial condition
T0      =  1300;                     % reference/top potential temperature [C]
T1      =  1350;                     % bottom potential temperature (if different from top) [C]

Ttype   = 'gaussian';                % set initial temperature field type

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
BCTTop      = 'insulating';         % 'isothermal' or 'insulating' bottom boundaries
BCTBot      = 'insulating';         % 'isothermal' or 'insulating' bottom boundaries
BCTSides    = 'insulating';         % 'isothermal' or 'insulating' bottom boundaries

% Velocity boundary conditions: free slip = -1; no slip = 1
BCsides     = -1;                   % side boundaries
BCtop       = -1;                   % top boundary
BCbot       = -1;                   % bottom boundary

%% set solver options
% advection scheme
ADVN        =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
BCA         = {'closed','periodic'};% boundary condition on advection (top/bot, sides)
TINY        = 1e-16;                % tiny number to safeguard [0,1] limits
lambda      = 1/2;   	            % iterative lagging for phase fractionCFL         = 0.25;   	            % Courant number to limit physical time step
reltol    	= 1e-6;                 % relative residual tolerance for nonlinear iterations
abstol      = 1e-9;                 % absolute residual tolerance for nonlinear iterations
maxit       = 30;                   % maximum iteration count
tauR        = 1e16;
CFL         = 1;                % (physical) time stepping courant number (multiplies stable step) [0,1]
etareg      = 1e0;                  % regularisation factor for viscosity
TINT        =  'bd3i';              % time integration scheme ('bwei','cnsi','bd3i','bd3s')

%% test time stepping
DDT = [h/2, h/4, h/8, h/16];

for dt = DDT
    %% set model timing
    yr              =  3600*24*365.25;       % seconds per year
    dtmax           =  dt;              % maximum time step
    maxstep         =  D/dt;                  % maximum number of time steps

    % [do not modify]
    tend            =  D/dt*h;           % model stopping time [s]
    % set temperature initial condition

rT      =  D/6;                  % radius of hot plume [m]
zT      =  D*0.5;                % z-position of hot plume [m]
xT      =  L/2;                  % x-position of hot plume [m]

    
    %% start model
    initialise_bnchm;
    Tin = T; rhoin = rho; Sin = S; XFein = XFe; XSiin = XSi; CFein = CFe; CSiin = CSi; xFein = xFe; xSiin = xSi;

    while time <= tend && step <= maxstep
        % print time step header
        fprintf(1,'\n*****  step = %d;  dt = %1.4e;  time = %1.4e yr;  %s\n\n',step,dt/yr,time/yr,dtlimit);

        figure(100);clf

    if     strcmp(TINT,'bwei') || step==1 % first step / 1st-order backward-Euler implicit scheme
        alpha1 = 1; alpha2 = 1; alpha3 = 0;
        beta1  = 1; beta2  = 0; beta3  = 0;
    elseif strcmp(TINT,'bd2i') || step==2 % second step / 2nd-order 3-point backward-difference implicit scheme
        alpha1 = 3/2; alpha2 = 4/2; alpha3 = -1/2;
        beta1  = 1;   beta2  =  0;  beta3  = 0;
    elseif strcmp(TINT,'cnsi')            % other steps / 2nd-order Crank-Nicolson semi-implicit scheme
        alpha1 = 1;   alpha2 = 1;   alpha3 = 0;
        beta1  = 1/2; beta2  = 1/2; beta3  = 0;
    elseif strcmp(TINT,'bd2s')            % other steps / 2nd-order 3-point backward-difference semi-implicit scheme
        alpha1 = 3/2; alpha2 = 4/2; alpha3 = -1/2;
        beta1  = 3/4; beta2  = 2/4; beta3  = -1/4;
    end

        % store previous solution and auxiliary fields
        Soo     = So;       So    = S;
        XFeoo   = XFeo;     XFeo  = XFe;
        XSioo   = XSio;     XSio  = XSi;
        CSioo   = CSio;     CSio  = CSi;
        CFeoo   = CFeo;     CFeo  = CFe;
        FsFeoo   = FsFeo;   FsFeo  = FsFe;
        FsSioo   = FsSio;   FsSio  = FsSi;
        dSdtoo  = dSdto;    dSdto = dSdt;
        dXFedtoo= dXFedto;  dXFedto = dXFedt;
        dXSidtoo= dXSidto;  dXSidto = dXSidt;
        dCFedtoo= dCFedto;  dCFedto = dCFedt;
        dCSidtoo= dCSidto;  dCSidto = dCSidt;
        dFFedtoo= dFFedto;  dFFedto = dFFedt;
        dFSidtoo= dFSidto;  dFSidto = dFSidt;
        rhooo   = rhoo;     rhoo    = rho;
        Div_rhoVoo = Div_rhoVo; Div_rhoVo = Div_rhoV;
        Div_Vo  = Div_V;
        dto     = dt;
        % temp
        cFeo = cFe; cSio = cSi;

        % temp for radioactive decay
        NAlo  = NAl;
        dNdto = dNdt;

        % reset residuals and iteration count
        resnorm   = 1e3;
        resnorm0  = resnorm;
        iter      = 0;


        % non-linear iteration loop
        while resnorm/resnorm0 >= reltol && resnorm >= abstol && iter <= maxit

            % solve thermo-chemical equations
            solve_thermochem;

            % update non-linear parameters and auxiliary variables
            up2date;
            W(:) = 0;               % z-velocity on z-face nodes
            U(:) = 1;               % x-velocity on x-face nodes
            P(:) = 0;  
            segsSi(:) = 0;
            segsFe(:) = 0;
            seglSi(:) = 0;
            seglFe(:) = 0;
            ks(:)     = 0; % update phase velocities
            WlSi = W + seglSi;
            UlSi = U;
            WsSi = W + segsSi;
            UsSi = U;
            WlFe = W + seglFe;
            UlFe = U;
            WsFe = W + segsFe;
            UsFe = U;


            if ~bnchm
                resnorm = resnorm_TC + resnorm_VP;
                if iter == 0
                    resnorm0 = resnorm+TINY;
                end
                fprintf(1,'  ---  it = %d;  abs res = %1.4e;  rel res = %1.4e  \n',iter,resnorm,resnorm/resnorm0)

                figure(100); if iter==1; clf; else; hold on; end
                plot(iter,log10(resnorm_TC),'b.',iter,log10(resnorm_VP),'r.',iter,log10(resnorm),'k.','MarkerSize',15,'LineWidth',1.5); box on; axis tight;
                drawnow;
            end

            iter = iter+1;
        end

        % update mass and energy errors
        history;


        fprintf(1,'         min T   =  %4.1f;    mean T   = %4.1f;    max T   = %4.1f;   [degC]\n' ,min(T(:)-273.15),mean(T(:)-273.15),max(T(:)-273.15));
        fprintf(1,'         min x_{Fe}   =  %1.4f;    mean x_{Fe}   = %1.4f;    max x_{Fe}   = %1.4f;   [wt]\n'   ,min(xFe(:)  ),mean(xFe(:)  ),max(xFe(:)  ));
        fprintf(1,'         min c_{Fe}   =  %1.4f;    mean c_{Fe}   = %1.4f;    max c_{Fe}   = %1.4f;   [wt]\n'   ,min(cFe(:)  ),mean(cFe(:)  ),max(cFe(:)  ));
        fprintf(1,'         min c_{Si}   =  %1.4f;    mean c_{Si}   = %1.4f;    max c_{Si}   = %1.4f;   [wt]\n'   ,min(cSi(:)  ),mean(cSi(:)  ),max(cSi(:)  ));

        fprintf(1,'         min sumX=  %4.1f;    mean sumX= %4.1f;    max sumX= %4.1f;   [kg/m3]\n'  ,min(XFe(:)+XSi(:)),mean(XFe(:)+XSi(:))   ,max(XFe(:)+XSi(:)));
        fprintf(1,'         min rho =  %4.1f;    mean rho = %4.1f;    max rho = %4.1f;   [kg/m3]\n'  ,min(rho(:)),mean(rho(:))   ,max(rho(:)));
        if ~mod(step,nop) %round(2*nop/CFL))
            output;
            %        suppfigs; % supplementary figures for testing
            if radheat
                figure(21); plot(HST.time/yr, HIST.Hr(2:end)); xlabel ('time'); ylabel('Hr (W/m^3)')
            end
        end

        % break script for NaN
        if isnan(resnorm)
            output
            print('Error, Breaking script')
            %         output
            return
        end
        % increment time
        dt;
        step = step + 1;
        time = time + dt;
        figure(101)
        plot(xP(2:end-1),Tin(N/2,2:end-1),'--k',xP(2:end-1),T(N/2,2:end-1),'-r'); axis ij tight; box on;
    end



    %% plot convergence
    EM = norm(rho(inz,inx)-rhoin(inz,inx),'fro')./norm(rhoin(inz,inx),'fro');
    ES = norm(S(inz,inx)-Sin(inz,inx),'fro')./norm(Sin(inz,inx),'fro');
    EXFe = norm(XFe(inz,inx)-XFein(inz,inx),'fro')./norm(XFein(inz,inx),'fro');
    EXSi = norm(XSi(inz,inx)-XSiin(inz,inx),'fro')./norm(XSiin(inz,inx),'fro');
    ECFe = norm(CFe(inz,inx)-CFein(inz,inx),'fro')./norm(CFein(inz,inx),'fro');
    ECSi = norm(CSi(inz,inx)-CSiin(inz,inx),'fro')./norm(CSiin(inz,inx),'fro');

    fh14 = figure(14);
    subplot(6,1,1);
    plot((maxstep/2+1/2)*dt:dt:(maxstep-1/2)*dt,diff(HST.EM(maxstep/2:maxstep))/dt,'-',LW{:}); hold on; axis tight; box on;
    ylabel('$\Delta E_M$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(6,1,2);
    plot((maxstep/2+1/2)*dt:dt:(maxstep-1/2)*dt,diff(HST.ES(maxstep/2:maxstep))/dt,'-',LW{:}); hold on; axis tight; box on;
    ylabel('$\Delta E_S$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(6,1,3);
    plot((maxstep/2+1/2)*dt:dt:(maxstep-1/2)*dt,diff(HST.EXFe(maxstep/2:maxstep))/dt,'-',LW{:}); hold on; axis tight; box on;
    ylabel('$\Delta E_XFe$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(6,1,4);
    plot((maxstep/2+1/2)*dt:dt:(maxstep-1/2)*dt,diff(HST.EXSi(maxstep/2:maxstep))/dt,'-',LW{:}); hold on; axis tight; box on;
    ylabel('$\Delta E_XSi$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(6,1,5);
    plot((maxstep/2+1/2)*dt:dt:(maxstep-1/2)*dt,diff(HST.ECFe(maxstep/2:maxstep))/dt,'-',LW{:}); hold on; axis tight; box on;
    ylabel('$\Delta E_CFe$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    xlabel('Time [s]',TX{:},FS{:});
    subplot(6,1,6);
    plot((maxstep/2+1/2)*dt:dt:(maxstep-1/2)*dt,diff(HST.ECSi(maxstep/2:maxstep))/dt,'-',LW{:}); hold on; axis tight; box on;
    ylabel('$\Delta E_CSi$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    xlabel('Time [s]',TX{:},FS{:});

    fh15 = figure(15);
    p1 = loglog(dt,EM,'kd','MarkerSize',8,'LineWidth',2); hold on; box on;
    p2 = loglog(dt,ES,'rs','MarkerSize',8,'LineWidth',2);
    p3 = loglog(dt,EXFe,'go','MarkerSize',8,'LineWidth',2);
    p4 = loglog(dt,EXSi,'bx','MarkerSize',8,'LineWidth',2);
    p5 = loglog(dt,ECFe,'m.','MarkerSize',8,'LineWidth',2);
    p6 = loglog(dt,ECSi,'c*','MarkerSize',8,'LineWidth',2);
    set(gca,'TicklabelInterpreter','latex','FontSize',12)
    xlabel('time stepping [dt]','Interpreter','latex','FontSize',16)
    ylabel('rel. numerical error rate [1/s]','Interpreter','latex','FontSize',16)
    title('Global conservation in time','Interpreter','latex','FontSize',20)

    if dt == DDT(1)
        p7 = loglog(DDT,ES.*(DDT./DDT(1)).^2,'k-','LineWidth',2);  % plot trend for comparison
    end
    if dt == DDT(end)
        legend([p1,p2,p3,p4,p5,p6,p7],{'error $M$','error $S$','error $X_{Fe}$','error $X_{Si}$','error $C_{Fe}$','error $C_{Si}$','2nd order'},'Interpreter','latex','box','on','location','southeast')
    end
    drawnow;
end

name = [outpath,'/',RunID,'_bnchm'];
print(fh15,name,'-dpng','-r300','-vector');