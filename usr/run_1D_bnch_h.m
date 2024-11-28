% planetesimal sill rainfall: user control script
% no sticky air/space; no self gravity
% equal grid spacing
clear ; close all
RunID           =  'bnchm_h_4phs_bd3i_dt_new3';               % run identifier
% create output directory
[~,systemname]  = system('hostname');
systemname(end) = [];
switch systemname
    case 'Horatio'
        outpath = ['/media/43TB_RAID_Array/fbintang/test_out/out/benchmark_h/',RunID];
        if ~exist(outpath, 'dir'); mkdir(outpath); end
    otherwise
        outpath = ['../out/benchmark_h/',RunID];
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
nop             =  10;                   % output every 'nop' grid steps of transport
bnchm           =  0;                    % manufactured solution benchmark on fluid mechanics solver
selfgrav = 0;
mode = 'cartesian';
%temporary
radheat         =  0;                    % radiogenic heating

D               =  500;                  % domain depth

%% set thermochemical parameters
% set initial system and component fractions
xFe0            =  0.5;                 % Fe-FeS system fraction
cFe0            =  0.15;                 % Fe-FeS fertile component fraction ([wt% S], maximum 0.35 for pure FeS
cSi0            =  0.42;                 % Si system fertile component fraction [wt% SiO2]

% set parameters
dxFe            =  0.1;                 % amplitude of initial random perturbation to iron system
dcFe            =  0.01;                 % amplitude of initial random perturbation to iron component
dcSi            =  0.01;                 % amplitude of initial random perturbation to silicate component

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
% rT      =  0;                  % radius of hot plume [m]
% zT      =  D*0.5;                % z-position of hot plume [m]
% xT      =  L/2;                  % x-position of hot plume [m]

Ttype   = 'gaussian';                % set initial temperature field type


%% set material parameters
% buoyancy parameters
rhosSi0     =  3300;                 % reference density solid refractory silicate [kg/m3]
rholSi0     =  2900;                 % reference density liquid refractory silicate [kg/m3]
rhosFe0     =  8000;                 % reference desnity solid refractory iron [kg/m3]
rholFe0     =  7600;                 % reference desnity liquid refractory iron [kg/m3]
gCSi        =  0.50;                 % compositional expansivity silicate
gCFe        =  0.65;                 % compositional expansivity iron
aT          =  0;                 % thermal expansivity silicate [1/K]
dx0         =  1e-3;                 % solid grain size [m]
df0         =  1e-3;                 % metal droplet size [m]
dm0         =  1e-3;                 % melt film size [m]
gz0         =  0.1;                  % z-gravity
gx0         =  0;               	 % x-gravity

% Reference pressure
P0 = 0;

% rheology parameters
EtalSi0     =  1e2;                  % reference silicate melt viscosity [Pas]
EtalFe0     =  1e1;                  % reference metal melt viscosity [Pas]
EtasSi0     =  1e15;                    % reference silicate/iron crystal viscosity
EtasFe0     =  1e19;
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


% segregation boundary conditions: 0 = depletion/accumulation; 1 =
% supply/sink
BCsegTop = 0;
BCsegBot = 0;

%% set solver options
% advection scheme
ADVN        =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
BCA         = {'periodic','closed'};% boundary condition on advection (top/bot, sides)
TINY        = 1e-16;                % tiny number to safeguard [0,1] limits
lambda      = 0.5;   	            % iterative lagging for phase fractionCFL         = 0.25;   	            % Courant number to limit physical time step
reltol    	= 1e-6;                 % relative residual tolerance for nonlinear iterations
abstol      = 1e-9;                 % absolute residual tolerance for nonlinear iterations
maxit       = 30;                   % maximum iteration count
alpha       = 0.50;                    % iterative step size parameter
beta        = 0.05;                    % iterative damping parameter
tauR        = 1e16;
CFL         = 1;                 % (physical) time stepping courant number (multiplies stable step) [0,1]
etareg      = 1e0;                  % regularisation factor for viscosity
TINT        =  'bd3i';              % time integration scheme ('bwei','cnsi','bd3i','bd3s')
mixReg = 0; SMALL = 1e-6;
dscale = 0; kmin = 0; etamin = EtalSi0;
%% test nonlinear tolerance
NN = [50, 100, 200];

for N = NN
    % [do not modify]
    h               =  D/N;          % spacing of x/z  coordinates
    L               =  h;
    
    % set gaussian positioning
    rT      =  D/6;                  % radius of hot plume [m]
    zT      =  D/2;                % z-position of hot plume [m]
    xT      =  L/2;                  % x-position of hot plume [m]
    %% set model timing
        % [do not modify]
    dt              =  D/NN(3)/16;           % (initial) time step [s]

    tend            =  D/dt*h;              % model stopping time [s]
    yr              =  3600*24*365.25;      % seconds per year
    dtmax           =  dt;                   % maximum time step
    maxstep         =  D/dt;                % maximum number of time steps


   

    %% start model
 initialise_bnchm;
    Tin = T; rhoin = rho; Sin = S; FlFein = FlFe; FsFein = FsFe; FlSiin = FlSi; FsSiin = FsSi; XSiin = XSi; CFein = CFe; CSiin = CSi; xFein = xFe; xSiin = xSi;

    while time <= tend && step <= maxstep
        % print time step header
        fprintf(1,'\n*****  step = %d;  dt = %1.4e;  time = %1.4e yr;  %s\n\n',step,dt/yr,time/yr,dtlimit);

        figure(100);clf

    if     strcmp(TINT,'bwei') || step==1 % first step / 1st-order backward-Euler implicit scheme
        a1 = 1; a2 = 1; a3 = 0;
        b1  = 1; b2  = 0; b3  = 0;
    elseif strcmp(TINT,'bd2i') || step==2 % second step / 2nd-order 3-point backward-difference implicit scheme
        a1 = 3/2; a2 = 4/2; a3 = -1/2;
        b1  = 1;   b2  =  0;  b3  = 0;
    elseif strcmp(TINT,'cnsi')            % other steps / 2nd-order Crank-Nicolson semi-implicit scheme
        a1 = 1;   a2 = 1;   a3 = 0;
        b1  = 1/2; b2  = 1/2; b3  = 0;
    elseif strcmp(TINT,'bd2s')            % other steps / 2nd-order 3-point backward-difference semi-implicit scheme
        a1 = 3/2; a2 = 4/2; a3 = -1/2;
        b1  = 3/4; b2  = 2/4; b3  = -1/4;
    end

        % store previous solution and auxiliary fields
        Soo         = So;       So          = S;
        XFeoo       = XFeo;     XFeo        = XFe;
        XSioo       = XSio;     XSio        = XSi;
        CSioo       = CSio;     CSio        = CSi;
        CFeoo       = CFeo;     CFeo        = CFe;
        FsFeoo      = FsFeo;    FsFeo       = FsFe;
        FsSioo      = FsSio;    FsSio       = FsSi;
        FlFeoo      = FlFeo;    FlFeo       = FlFe;
        FlSioo      = FlSio;    FlSio       = FlSi;
        dSdtoo      = dSdto;    dSdto       = dSdt;
        dXFedtoo    = dXFedto;  dXFedto     = dXFedt;
        dXSidtoo    = dXSidto;  dXSidto     = dXSidt;
        dCFedtoo    = dCFedto;  dCFedto     = dCFedt;
        dCSidtoo    = dCSidto;  dCSidto     = dCSidt;
        dFsFedtoo   = dFsFedto; dFsFedto    = dFsFedt;
        dFsSidtoo   = dFsSidto; dFsSidto    = dFsSidt;
        dFlFedtoo   = dFlFedto; dFlFedto    = dFlFedt;
        dFlSidtoo   = dFlSidto; dFlSidto    = dFlSidt;
        rhooo       = rhoo;     rhoo        = rho;
        advn_RHOoo  = advn_RHOo;advn_RHOo   = advn_RHO;
        Div_Vo      = Div_V;
        dto         = dt;
    To = T;
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
            % solve_thermochem;

            % store previous iteration (diagnostic)
            Ti    = T;
            Si    = S;
            XFei  = XFe;
            XSii  = XSi;
            CFei  = CFe;
            CSii  = CSi;
            FlFei = FlFe;
            FlSii = FlSi;
            %% update heat content (entropy)
            advn_S = - advect(rp(inz,:).^2.*FsSi(inz,inx).*ssSi(inz,inx),UsSi(inz,:),WsSi(:,inx),h,{ADVN,''},[1,2],BCA)./rp(inz,:).^2 ...  % heat advection
                - advect(rp(inz,:).^2.*FsFe(inz,inx).*ssFe(inz,inx),UsFe(inz,:),WsFe(:,inx),h,{ADVN,''},[1,2],BCA)./rp(inz,:).^2 ...
                - advect(rp(inz,:).^2.*FlSi(inz,inx).*slSi(inz,inx),UlSi(inz,:),WlSi(:,inx),h,{ADVN,''},[1,2],BCA)./rp(inz,:).^2 ...
                - advect(rp(inz,:).^2.*FlFe(inz,inx).*slFe(inz,inx),UlFe(inz,:),WlFe(:,inx),h,{ADVN,''},[1,2],BCA)./rp(inz,:).^2;
            dSdt   = advn_S;
            % residual of entropy evolution
            res_S = (a1*S(inz,inx)-a2*So(inz,inx)-a3*Soo(inz,inx))/dt - (b1*dSdt + b2*dSdto + b3*dSdtoo);
            % update solution
            upd_S       = - alpha*res_S*dt/a1 + beta*upd_S;
            S(inz,inx)  = S(inz,inx) + upd_S;
            S(1,:)  = S(2,:); S(end,:) = S(end-1,:); S(:,[1 end]) = S(:,[2 end-1]);

            %% update composition

            % update fertile chemical components
            advn_CSi = - advect(rp(inz,:).^2.*FsSi(inz,inx).*csSi(inz,inx),UsSi(inz,:),WsSi(:,inx),h,{ADVN,''},[1,2],BCA)./rp(inz,:).^2 ...
                       - advect(rp(inz,:).^2.*FlSi(inz,inx).*clSi(inz,inx),UlSi(inz,:),WlSi(:,inx),h,{ADVN,''},[1,2],BCA)./rp(inz,:).^2;
            dCSidt   = advn_CSi;

            advn_CFe = - advect(rp(inz,:).^2.*FsFe(inz,inx).*csFe(inz,inx),UsFe(inz,:),WsFe(:,inx),h,{ADVN,''},[1,2],BCA)./rp(inz,:).^2 ...
                       - advect(rp(inz,:).^2.*FlFe(inz,inx).*clFe(inz,inx),UlFe(inz,:),WlFe(:,inx),h,{ADVN,''},[1,2],BCA)./rp(inz,:).^2;
            dCFedt   = advn_CFe + diff_CFe;

            % update solution
            % residual of entropy evolution
            res_CFe = (a1*CFe(inz,inx)-a2*CFeo(inz,inx)-a3*CFeoo(inz,inx))/dt - (b1*dCFedt + b2*dCFedto + b3*dCFedtoo);
            res_CSi = (a1*CSi(inz,inx)-a2*CSio(inz,inx)-a3*CSioo(inz,inx))/dt - (b1*dCSidt + b2*dCSidto + b3*dCSidtoo);

            % update solution
            % semi-implicit update of bulk chemical composition density
            upd_CFe       = - alpha*res_CFe*dt/a1 + beta*upd_CFe;
            upd_CSi       = - alpha*res_CSi*dt/a1 + beta*upd_CSi;
            CFe(inz,inx)  = CFe(inz,inx) + upd_CFe; CSi(inz,inx)  = CSi(inz,inx) + upd_CSi;

            % apply boundaries
            CSi([1 end],:) = CSi([2 end-1],:);  CSi(:,[1 end]) = CSi(:,[2 end-1]);
            CFe([1 end],:) = CFe([2 end-1],:);  CFe(:,[1 end]) = CFe(:,[2 end-1]);

            % % enforce 0,rho limits
            CFe = max(0,CFe); CSi = max(0,CSi);

            %% update local phase equilibrium
            [fsFeq,csFeq,clFeq] = equilibrium(T,cFe,Pt,TFe1,TFe2,cphsFe1,cphsFe2,...
                perTFe,percsFe,perclFe,clap,PhDgFe);
            % apply boundaries
            fsFeq([1 end],:) = fsFeq([2 end-1],:); fsFeq(:,[1 end]) = fsFeq(:,[2 end-1]);
            csFeq([1 end],:) = csFeq([2 end-1],:); csFeq(:,[1 end]) = csFeq(:,[2 end-1]);
            clFeq([1 end],:) = clFeq([2 end-1],:); clFeq(:,[1 end]) = clFeq(:,[2 end-1]);
            [fsSiq,csSiq,clSiq] = equilibrium(T,cSi,Pt,TSi1,TSi2,cphsSi1,cphsSi2,...
                perTSi,percsSi,perclSi,clap,PhDgSi);
            % apply boundaries
            fsSiq([1 end],:) = fsSiq([2 end-1],:); fsSiq(:,[1 end]) = fsSiq(:,[2 end-1]);
            csSiq([1 end],:) = csSiq([2 end-1],:); csSiq(:,[1 end]) = csSiq(:,[2 end-1]);
            clSiq([1 end],:) = clSiq([2 end-1],:); clSiq(:,[1 end]) = clSiq(:,[2 end-1]);
            flFeq = 1-fsFeq; flSiq = 1-fsSiq;
            
            %% update phase fractions
            % solid
            GsFe        = ((XFe.*fsFeq-FsFe)./(4.*dt));
            GlFe        = ((XFe.*flFeq-FlFe)./(4.*dt));
            GsSi        = ((XSi.*fsSiq-FsSi)./(4.*dt));
            GlSi        = ((XSi.*flSiq-FlSi)./(4.*dt));
            advn_FsFe   = - advect(rp(inz,:).^2.*FsFe(inz,inx),UsFe(inz,:),WsFe(:,inx),h,{ADVN,''},[1,2],BCA)./rp(inz,:).^2;
            dFsFedt     = advn_FsFe + GsFe(inz,inx);
            res_FsFe    = (a1*FsFe(inz,inx)-a2*FsFeo(inz,inx)-a3*FsFeoo(inz,inx))/dt - (b1*dFsFedt + b2*dFsFedto + b3*dFsFedtoo);

            advn_FsSi   = - advect(rp(inz,:).^2.*FsSi(inz,inx),UsSi(inz,:),WsSi(:,inx),h,{ADVN,''},[1,2],BCA)./rp(inz,:).^2;
            dFsSidt     = advn_FsSi + GsSi(inz,inx);                                       % total rate of change
            res_FsSi    = (a1*FsSi(inz,inx)-a2*FsSio(inz,inx)-a3*FsSioo(inz,inx))/dt - (b1*dFsSidt + b2*dFsSidto + b3*dFsSidtoo);

            upd_FsFe      =   - alpha*res_FsFe*dt/a1 + beta*upd_FsFe;
            upd_FsSi      =   - alpha*res_FsSi*dt/a1 + beta*upd_FsSi;
            FsFe(inz,inx) = FsFe(inz,inx) + upd_FsFe;
            FsSi(inz,inx) = FsSi(inz,inx) + upd_FsSi;


            FsFe([1 end],:) = FsFe([2 end-1],:);  FsFe(:,[1 end]) = FsFe(:,[2 end-1]);  % apply boundary conditions
            FsSi([1 end],:) = FsSi([2 end-1],:);  FsSi(:,[1 end]) = FsSi(:,[2 end-1]);  % apply boundary conditions
            FsFe = max(0, FsFe ); FsSi = max(0, FsSi );
            

            % liquid
            advn_FlFe   = - advect(rp(inz,:).^2.*FlFe(inz,inx),UlFe(inz,:),WlFe(:,inx),h,{ADVN,''},[1,2],BCA)./rp(inz,:).^2;
            dFlFedt       = advn_FlFe + GlFe(inz,inx);
            res_FlFe = (a1*FlFe(inz,inx)-a2*FlFeo(inz,inx)-a3*FlFeoo(inz,inx))/dt - (b1*dFlFedt + b2*dFlFedto + b3*dFlFedtoo);

            advn_FlSi   = - advect(rp(inz,:).^2.*FlSi(inz,inx),UlSi(inz,:),WlSi(:,inx),h,{ADVN,''},[1,2],BCA)./rp(inz,:).^2;
            dFlSidt     = advn_FlSi + GlSi(inz,inx);                                       % total rate of change
            res_FlSi = (a1*FlSi(inz,inx)-a2*FlSio(inz,inx)-a3*FlSioo(inz,inx))/dt - (b1*dFlSidt + b2*dFlSidto + b3*dFlSidtoo);

            upd_FlFe       =   - alpha*res_FlFe*dt/a1 + beta*upd_FlFe;
            upd_FlSi       =   - alpha*res_FlSi*dt/a1 + beta*upd_FlSi;
            FlFe(inz,inx)  = FlFe(inz,inx) + upd_FlFe;
            FlSi(inz,inx)  = FlSi(inz,inx) + upd_FlSi;

            FlFe([1 end],:) = FlFe([2 end-1],:);  FlFe(:,[1 end]) = FlFe(:,[2 end-1]);  % apply boundary conditions
            FlSi([1 end],:) = FlSi([2 end-1],:);  FlSi(:,[1 end]) = FlSi(:,[2 end-1]);  % apply boundary conditions

            FlFe = max(0, FlFe ); FlSi = max(0, FlSi );            
            XFe  = FsFe + FlFe; XSi  = FsSi + FlSi;
            
            RHO  = XFe + XSi;
            advn_RHO = advn_FsFe + advn_FlFe + advn_FsSi + advn_FlSi;

            % % update temperature
            T   = T0.*exp((S - FsFe.*dEntrFe - FsSi.*dEntrSi)./(FlFe+FsFe+FlSi+FsSi)./Cp ...
                + aT.*(Pt - P0)./rhoRef./Cp);
            Tp  = T0.*exp((S - FsFe.*dEntrFe - FsSi.*dEntrSi)./(FlFe+FsFe+FlSi+FsSi)./Cp);

            % update system fractions
            xFe = max(0,min(1, XFe./RHO));
            xSi = max(0,min(1, XSi./RHO ));

            hasFe   = xFe>SMALL;
            hasSi   = xSi>SMALL;
            if any(hasFe<=SMALL)
                hasFe;
            end

            % update chemical composition
            cSi(hasSi) = CSi(hasSi)./XSi(hasSi);
            cFe(hasFe) = CFe(hasFe)./XFe(hasFe);
            cSi(~hasSi) = fsSiq(~hasSi).*csSiq(~hasSi) + flSiq(~hasSi).*clSiq(~hasSi);
            cFe(~hasFe) = fsFeq(~hasFe).*csFeq(~hasFe) + flFeq(~hasFe).*clFeq(~hasFe);

            % ensure real numbers
            cSi(isnan(cSi)) = 0; cFe(isnan(cFe)) = 0;

            % update phase fractions [wt]
            flFe(hasFe) = max(0,min(1,FlFe(hasFe)./XFe(hasFe)));
            flSi(hasSi) = max(0,min(1,FlSi(hasSi)./XSi(hasSi)));
            fsFe(hasFe) = max(0,min(1,FsFe(hasFe)./XFe(hasFe)));
            fsSi(hasSi) = max(0,min(1,FsSi(hasSi)./XSi(hasSi)));
            % ensure real numbers
            flFe(isnan(flFe)) = 0; fsFe(isnan(fsFe)) = 0;
            flSi(isnan(flSi)) = 0; fsFe(isnan(fsSi)) = 0;

            % flFe(hasFe) = 1-fsFe(hasFe);
            % flSi(hasSi) = 1-fsSi(hasSi);

            %detect where is fully solid or molten
            hassolSi = fsSi>SMALL;
            hassolFe = fsFe>SMALL;

            hasliqSi = flSi>SMALL;
            hasliqFe = flFe>SMALL;

            %% update phase compositions
            KcFe = csFeq./clFeq;
            clFe(hasFe) = max(TINY,cFe(hasFe))./max(SMALL,(flFe(hasFe) + fsFe(hasFe).*KcFe(hasFe)));
            csFe(hasFe) = max(TINY,cFe(hasFe))./max(SMALL,(flFe(hasFe)./KcFe(hasFe) + fsFe(hasFe)));
            ind = ~hassolFe | ~hasliqFe;
            csFe(ind) = csFeq(ind);
            clFe(ind) = clFeq(ind);

            KcSi = csSiq./clSiq;
            clSi(hasSi) = max(TINY,cSi(hasSi)./max(SMALL,(flSi(hasSi) + fsSi(hasSi).*KcSi(hasSi))));
            csSi(hasSi) = max(TINY,cSi(hasSi)./max(SMALL,(flSi(hasSi)./KcSi(hasSi) + fsSi(hasSi))));
            ind = ~hassolSi | ~hasliqSi;
            csSi(ind) = csSiq(ind);
            clSi(ind) = clSiq(ind);

            %% update phase entropies
            slFe  = (S - FsFe.*dEntrFe - FsSi.*dEntrSi)./RHO;
            slSi  = slFe;
            ssFe  = slFe + dEntrFe;
            ssSi  = slSi + dEntrSi;

            %% get residual of thermochemical equations from iterative update
            normS   = norm(S    - Si   ,2)./(norm(S   ,2)+TINY);
            normXFe = norm(XFe  - XFei ,2)./(norm(XFe ,2)+TINY);
            normXSi = norm(XSi  - XSii ,2)./(norm(XSi ,2)+TINY);
            normCFe = norm(CFe  - CFei ,2)./(norm(CFe ,2)+TINY);
            normCSi = norm(CSi  - CSii ,2)./(norm(CSi ,2)+TINY);
            normFFe = norm(FlFe - FlFei,2)./(norm(FlFe,2)+TINY);
            normFSi = norm(FlSi - FlSii,2)./(norm(FlSi,2)+TINY);
            resnorm_TC = normS + normXFe + normXSi + normCSi + normCFe + normFFe + normFSi;


            % update non-linear parameters and auxiliary variables
            % up2date;
            %% update T and C-dependent density
            rhosSi = rhosSi0.*(1 - aT.*(T-perTSi) - gCSi.*(csSi-cphsSi1));
            rholSi = rholSi0.*(1 - aT.*(T-perTSi) - gCSi.*(clSi-cphsSi1));
            rhosFe = rhosFe0.*(1 - aT.*(T-perTFe) - gCFe.*(csFe-cphsFe1));
            rholFe = rholFe0.*(1 - aT.*(T-perTFe) - gCFe.*(clFe-cphsFe1));

            % update mixture density
            rho    = 1./(xFe.*(fsFe./rhosFe + flFe./rholFe) ...
                + xSi.*(fsSi./rhosSi + flSi./rholSi));
            rho([1 end],:) = rho([2 end-1],:);  rho(:,[1 end]) = rho(:,[2 end-1]);

            phisFe = max(0,min(1,xFe.* fsFe .* rho ./ rhosFe));
            philFe = max(0,min(1,xFe.* flFe .* rho ./ rholFe));
            phisSi = max(0,min(1,xSi.* fsSi .* rho ./ rhosSi));
            philSi = max(0,min(1,xSi.* flSi .* rho ./ rholSi));

            
            W(:) = 1;               % z-velocity on z-face nodes
            U(:) = 0;               % x-velocity on x-face nodes
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
            TX = {'Interpreter','Latex'}; FS = {'FontSize',12};
            TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',10};
            UN = {'Units','Centimeters'};
            CL = {'Color',[0.0 0.0 0.0],[0.80 0.15 0.10],[0.10 0.15 0.65],[0.45 0.60 0.95]};
            LW = {'LineWidth',2};
            fh1 = figure(1); clf;
            sgtitle(['time = ',num2str(time/yr,3),' [yr]'],TX{:},FS{:},'Color','k');
            subplot(1,4,1)
            plot(mean(Tp(2:end-1,2:end-1),2)-273.15,zP(2:end-1)./1000,'LineWidth',2); axis ij tight; box on;
            title('$T [^\circ$C]',TX{:},FS{:}); set(gca,TL{:},TS{:});
            ylabel('Depth [km]',TX{:},FS{:});
            subplot(1,4,2)
            plot(mean(xFe(2:end-1,2:end-1),2).*100,zP(2:end-1)./1000,'r-','LineWidth',2); hold on;  axis ij tight; box on;
            plot(mean(xSi(2:end-1,2:end-1),2).*100,zP(2:end-1)./1000,'b-','LineWidth',2);
            title('$x_{Fe}$ / $x_{Si}$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
            legend('$x_{\mathrm{Fe}}$', '$x_{\mathrm{Si}}$',TX{:},'location','west')
            subplot(1,4,3)
            plot(mean(clFe(2:end-1,2:end-1),2).*100,zP(2:end-1)./1000,'-r','LineWidth',2); axis ij tight; box on; hold on
            plot(mean(csFe(2:end-1,2:end-1),2).*100,zP(2:end-1)./1000,'-b','LineWidth',2);
            plot(mean(cFe(2:end-1,2:end-1),2).*100,zP(2:end-1)./1000,'k-','LineWidth',2);
            title('$\bar{c}_{Fe}$ [wt\% S]',TX{:},FS{:}); set(gca,TL{:},TS{:});
            subplot(1,4,4)
            plot((mean(clSi(2:end-1,2:end-1),2)+cSimin).*100,zP(2:end-1)./1000,'-r','LineWidth',2); axis ij tight; box on; hold on
            plot((mean(csSi(2:end-1,2:end-1),2)+cSimin).*100,zP(2:end-1)./1000,'-b','LineWidth',2);
            plot((mean(cSi(2:end-1,2:end-1),2)+cSimin).*100,zP(2:end-1)./1000,'k-','LineWidth',2);
            title('$\bar{c}_{Si}$ [wt\% SiO$_2$]',TX{:},FS{:}); set(gca,TL{:},TS{:});
            legend('melt', 'solid', 'bulk', TX{:},'location','west')

            fh2 = figure(2); clf;
            sgtitle(['time = ',num2str(time/yr,3),' [yr]'],TX{:},FS{:},'Color','k');
            subplot(1,4,1)
            plot(mean(flFe(2:end-1,2:end-1),2).*100,zP(2:end-1)./1000,'LineWidth',2); axis ij tight; box on; hold on
            plot(mean(flSi(2:end-1,2:end-1),2).*100,zP(2:end-1)./1000,'LineWidth',2);
            title('$f_{j}^\ell$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
            legend('Fe melt', 'Si melt', TX{:},'location','west')
            ylabel('Depth [km]',TX{:},FS{:});
            subplot(1,4,2)
            plot(mean(phisFe(2:end-1,2:end-1),2).*100,zP(2:end-1)./1000,'LineWidth',2); axis ij tight; box on; hold on
            plot(mean(philFe(2:end-1,2:end-1),2).*100,zP(2:end-1)./1000,'LineWidth',2);
            plot(mean(phisSi(2:end-1,2:end-1),2).*100,zP(2:end-1)./1000,'LineWidth',2);
            plot(mean(philSi(2:end-1,2:end-1),2).*100,zP(2:end-1)./1000,'LineWidth',2);
            title('$\phi_{j}^i$ [vol\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
            legend('$\phi_{Fe}^s$','$\phi_{Fe}^{l}$','$\phi_{Si}^s$','$\phi_{Si}^{l}$',TX{:},'location','west')
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
        plot(mean(Tin(2:end-1,2:end-1),2),zP(2:end-1),'--k',mean(T(2:end-1,2:end-1),2),zP(2:end-1),'-r'); axis ij tight; box on;

    end



    %% plot convergence
    EM = norm(rho(inz,inx)-rhoin(inz,inx),'fro')./norm(rhoin(inz,inx),'fro');
    ES = norm(S(inz,inx)-Sin(inz,inx),'fro')./norm(Sin(inz,inx),'fro');
    EFlFe = norm(FlFe(inz,inx)-FlFein(inz,inx),'fro')./norm(FlFein(inz,inx),'fro');
    EFsFe = norm(FsFe(inz,inx)-FsFein(inz,inx),'fro')./norm(FsFein(inz,inx),'fro');
    EFlSi = norm(FlSi(inz,inx)-FlSiin(inz,inx),'fro')./norm(FlSiin(inz,inx),'fro');
    EFsSi = norm(FsSi(inz,inx)-FsSiin(inz,inx),'fro')./norm(FsSiin(inz,inx),'fro');
    ECFe = norm(CFe(inz,inx)-CFein(inz,inx),'fro')./norm(CFein(inz,inx),'fro');
    ECSi = norm(CSi(inz,inx)-CSiin(inz,inx),'fro')./norm(CSiin(inz,inx),'fro');

    % fh14 = figure(14);
    % subplot(6,1,1);
    % plot((maxstep/2+1/2)*dt:dt:(maxstep-1/2)*dt,diff(HST.EM(maxstep/2:maxstep))/dt,'-',LW{:}); hold on; axis tight; box on;
    % ylabel('$\Delta E_M$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    % subplot(6,1,2);
    % plot((maxstep/2+1/2)*dt:dt:(maxstep-1/2)*dt,diff(HST.ES(maxstep/2:maxstep))/dt,'-',LW{:}); hold on; axis tight; box on;
    % ylabel('$\Delta E_S$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    % subplot(6,1,3);
    % plot((maxstep/2+1/2)*dt:dt:(maxstep-1/2)*dt,diff(HST.EXFe(maxstep/2:maxstep))/dt,'-',LW{:}); hold on; axis tight; box on;
    % ylabel('$\Delta E_XFe$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    % subplot(6,1,4);
    % plot((maxstep/2+1/2)*dt:dt:(maxstep-1/2)*dt,diff(HST.EXSi(maxstep/2:maxstep))/dt,'-',LW{:}); hold on; axis tight; box on;
    % ylabel('$\Delta E_XSi$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    % subplot(6,1,5);
    % plot((maxstep/2+1/2)*dt:dt:(maxstep-1/2)*dt,diff(HST.ECFe(maxstep/2:maxstep))/dt,'-',LW{:}); hold on; axis tight; box on;
    % ylabel('$\Delta E_CFe$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    % xlabel('Time [s]',TX{:},FS{:});
    % subplot(6,1,6);
    % plot((maxstep/2+1/2)*dt:dt:(maxstep-1/2)*dt,diff(HST.ECSi(maxstep/2:maxstep))/dt,'-',LW{:}); hold on; axis tight; box on;
    % ylabel('$\Delta E_CSi$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    % xlabel('Time [s]',TX{:},FS{:});

    fh15 = figure(15);
    p1 = loglog(h,EM,'kd','MarkerSize',8,'LineWidth',2); hold on; box on;
    p2 = loglog(h,ES,'rs','MarkerSize',8,'LineWidth',2);
    p3 = loglog(h,EFlFe,'go','MarkerSize',8,'LineWidth',2);
    p4 = loglog(h,EFsFe,'bx','MarkerSize',8,'LineWidth',2);
    p5 = loglog(h,EFlSi,'g+','MarkerSize',8,'LineWidth',2);
    p6 = loglog(h,EFsSi,'b*','MarkerSize',8,'LineWidth',2);
    p7 = loglog(h,ECFe,'m.','MarkerSize',8,'LineWidth',2);
    p8 = loglog(h,ECSi,'c*','MarkerSize',8,'LineWidth',2);
    set(gca,'TicklabelInterpreter','latex','FontSize',12)
    xlabel('grid spacing [m]','Interpreter','latex','FontSize',16)
    ylabel('rel. numerical error rate [1/s]','Interpreter','latex','FontSize',16)
    title('Global conservation in space','Interpreter','latex','FontSize',20)

    if N == NN(1)
        p9 = loglog(D./NN,ES.*((D./NN)./(D./NN(1))).^1,'k-','LineWidth',2); hold on  % plot trend for comparison
        p10 = loglog(D./NN,ES.*((D./NN)./(D./NN(1))).^2,'k-','LineWidth',2);
    end
    if N == NN(end)
        legend([p1,p2,p3,p4,p5,p6,p7,p8,p9,p10],{'error $M$','error $S$','error $F^l_{\mathrm{Fe}}$','error $F^s_{\mathrm{Fe}}$',...
            'error $F^l_{\mathrm{Si}}$','error $F^s_{\mathrm{Si}}$','error $C_{\mathrm{Fe}}$','error $C_{\mathrm{Si}}$','2nd order','3rd order'},...
            'Interpreter','latex','box','on','location','southeast')
    end
    drawnow;
end

name = [outpath,'/',RunID,'_bnchm'];
print(fh15,name,'-dpng','-r300','-vector');
savefig(name)