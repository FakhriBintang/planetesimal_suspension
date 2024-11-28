% suspension flow: initialise model run

% print initialisation header
fprintf(1,'  ---  initialise model run \n\n');

if save_op
copyfile(infile, outpath);
end

%% setup numerical grid
% set dimensions of staggered/ghosted 2D grid

nxC         =  Nx+1;                            % number of corner nodes in x-direction
nzC         =  Nz+1;                            % number of corner nodes in z-direction
nxP         =  Nx+2;                            % number of centre nodes in x-direction
nzP         =  Nz+2;                            % number of centre nodes in z-direction
nxW         =  Nx+2;                            % number of z-face nodes in x-direction
nzW         =  Nz+1;                            % number of z-face nodes in z-direction
nxU         =  Nx+1;                            % number of x-face nodes in x-direction
nzU         =  Nz+2;                            % number of x-face nodes in z-direction

switch mode
    case 'spherical'
        % change dimensionz to z by 1 for 1D spherical fluid mechanics
        NC = nzC; NP = nzP; NW = nzW; NU = nxU; NDOF = NP+NW+NU;
    otherwise
        NC          =  nxC*nzC;                         % total number of corner nodes
        NP          =  nxP*nzP;                         % total number of corner nodes
        NW          =  nxW*nzW;                         % total number of z-face nodes
        NU          =  nxU*nzU;                         % total number of x-face nodes
        NDOF        =  NP+NW+NU;                        % total number of all degrees of freedum
end

% set coordinate vectors
xC          =  0:h:L;                           % Horizontal coordinates of corner nodes [m]
zC          =  0:h:D;                           % Vertical   coordinates of corner nodes [m]
xP          =  -h/2:h:L+h/2;                    % Horizontal coordinates of centre nodes [m]
zP          =  -h/2:h:D+h/2;                    % Vertical   coordinates of centre nodes [m]
xW          =  xP;                              % Horizontal coordinates of z-face nodes [m]
zW          =  zC;                              % Vertical   coordinates of z-face nodes [m]
xU          =  xC;                              % Horizontal coordinates of x-face nodes [m]
zU          =  zP;                              % Vertical   coordinates of x-face nodes [m]

% set 2D coordinate grids
[XC,ZC] = meshgrid(xC,zC);              % corner nodes grid
[XP,ZP] = meshgrid(xP,zP);              % centre nodes grid
[XW,ZW] = meshgrid(xW,zW);              % z-face nodes grid
[XU,ZU] = meshgrid(xU,zU);              % x-face nodes grid


%% setup mapping arrays
%% apply radius if using spherical coordiantes
rp          =  ones(size(zP))';
rw          =  ones(size(zW))';
switch mode
    case'spherical'
        rp          =  flip(rmin-h/2:h:D+h/2+rmin)';
        rw          =  (rp(1:end-1)+rp(2:end))./2;
        % change dimensionz to z by 1 for 1D spherical fluid mechanics
        NC = nzC; NP = nzP; NW = nzW; NU = nxU; NDOF = NP+NW+NU;
        Map         =  reshape(1:NP,nzP,1);
        MapW        =  reshape(1:NW,nzW,1);
        % MapU        =  reshape(1:NU,nzU,1) + NW;
        MapP        =  reshape(1:NP,nzP,1);
    case'cartesian'
        Map         =  reshape(1:NP,nzP,nxP);
        MapW        =  reshape(1:NW,nzW,nxW);
        MapU        =  reshape(1:NU,nzU,nxU) + NW;
        MapP        =  reshape(1:NP,nzP,nxP);
end

inz = 2:length(zP)-1;
inx = 2:length(xP)-1;

% get smoothed initialisation field
rng(67);
pert = randn(nzP,nxP);
for i = 1:round(smth)
    pert(2:end-1,2:end-1) = pert(2:end-1,2:end-1) + diff(pert(:,2:end-1),2,1)./8 ...
                                              + diff(pert(2:end-1,:),2,2)./8;
    pert([1 end],:)       = pert([2 end-1],:);
    pert(:,[1 end])       = pert(:,[2 end-1]);
end
pert              = pert./max(abs(pert(:)));
pert              = pert - mean(mean(pert(2:end-1,2:end-1)));


%% setup material property arrays
% gravity
gxP = zeros(nzP,nxP) + gx0;      gzP = zeros(nzP,nxP) + gz0;
gx  = zeros(nzU,nxU) + gx0;      gz  = zeros(nzW,nxW) + gz0;


%% setup velocity-pressure solution arrays
W           = zeros(nzW,nxW);               % z-velocity on z-face nodes
WBG = W; %wf = 0.*W; wx = 0.*W; wm = 0.*W; upd_W = 0*W;
U           = zeros(nzU,nxU);               % x-velocity on x-face nodes
UBG = U;
P           = zeros(nzP,nxP);               % pressure on centre nodes
SOL         = [W(:);U(:);P(:)];             % full solution vector
Vel         = zeros(nzP,nxP);               % velocity magnitude

% projected velocities on centre nodes
UP              = zeros(nzP,nxP);
WP              = zeros(nzP,nxP);
UP(:,2:end-1)   = (U(:,1:end-1)+U(:,2:end))./2;
WP(2:end-1,:)   = (W(1:end-1,:)+W(2:end,:))./2;

%% setup deformation property arrays
ups = zeros(nzP,nxP);               % velocity divergence on centre nodes
exx = zeros(nzP,nxP);               % x-normal strain rate on centre nodes
ezz = zeros(nzP,nxP);               % z-normal strain rate on centre nodes
exz = zeros(nzC,nxC);               % xz-shear strain rate on corner nodes
eII = zeros(nzP,nxP);               % strain rate magnitude on centre nodes
txx = zeros(nzP,nxP);               % x-normal stress on centre nodes
tzz = zeros(nzP,nxP);               % z-normal stress on centre nodes
txz = zeros(nzC,nxC);               % xz-shear stress on corner nodes
tII = zeros(nzP,nxP);               % stress magnitude on centre nodes


%% setup heating rates
dndt        = 0;
dSdt        = zeros(Nz,Nx);             % entropy rate of change
dXFedt      = zeros(Nz,Nx);             % Iron system rate of change 
dXSidt      = zeros(Nz,Nx);             % Si system rate of change 
dCSidt      = zeros(Nz,Nx);             % Silicate component density rate of change
dCFedt      = zeros(Nz,Nx);             % Iron component density rate of change
dFsFedt     = zeros(Nz,Nx);             % Iron solid fraction rate of change
dFsSidt     = zeros(Nz,Nx);             % Silicate solid fraction rate of change
dFlFedt     = zeros(Nz,Nx);             % Iron liquid fraction rate of change
dFlSidt     = zeros(Nz,Nx);             % Silicate liquid fraction rate of change
diff_S      = zeros(Nz,Nx);             % Heat dissipation rate
diss_T      = zeros(Nz,Nx);             % Temperature diffusion rate
diff_CSi    = zeros(Nz,Nx);             % Silicate component diffusion rate
diff_CFe    = zeros(Nz,Nx);             % Iron component diffusion rate
Div_V       = zeros(Nz+2,Nx+2);         % Stokes velocity divergence
advn_RHO    = zeros(Nz,Nx);             % Mixture mass flux divergence
advn_FsFe   = zeros(Nz,Nx);             % iron solid fraction advection 
advn_FlFe   = zeros(Nz,Nx);             % iron liquid fraction advection 
advn_FsSi   = zeros(Nz,Nx);             % silicate solid fraction advection 
advn_FlSi   = zeros(Nz,Nx);             % silicate liquid fraction advection 
VolSrc      = zeros(Nz,Nx);             % volume source term

% initialise counting variables
RUN.frame   = 0;      % initialise output frame count
time        = 0;      % initialise time count
step        = 0;      % initialise time step count
iter        = 0;      % initialise iteration count

%% set initial condition on thermochemical fields
% temperature
pert = -h/2.*cos(XP*2*pi/D);
% temporary, set theral distribution parameters. Switch to run scripts
% later
zlay     =  0.0;                 % layer thickness (relative to domain depth D)
wlay_T   =  0.1;                % thickness of smooth layer boundary (relative to domain depth D)
% wlay_c   =  2*h/D;       % thickness of smooth layer boundary (relative to domain depth D)

switch Ttype
    case 'constant'     % constant temperature
        T      = T0 + (T1-T0) .* erf((ZP/D)/wlay_T); % + dT.*rp;
        T(1,:) = T0;
    case 'linear'       % linear temperature gradient with depth
        T      = T0 + abs(ZP+pert)./D.*(T1-T0);
    case 'gaussian'     % constant temperature with gaussian plume
        T      = zeros(nzP,nxP) + T0;

        T      = T + (T1-T0).*exp(-(XP-xT).^2./rT.^2 - (ZP-zT).^2./rT.^2 );
    case 'hot bottom'
        T      = zeros(nzP,nxP) + T0;
        T(end-10:end,:) = T0+100;
end
Tp = T;  % initial condition sets potential temperature [k]

% set initial component weight fraction [kg/kg]
xFe = xFe0 + dxFe.*pert; % Fe system
xSi = 1 - xFe;         % Si system
cFe = zeros(size(xFe)) + cFe0 + dcFe.*pert; % Fe component
cSi = zeros(size(xSi)) + cSi0 + dcSi.*pert - cSimin; % Si component
Pt = 0;

Hr  = Hr0.*zeros(size(T));

%% initialise loop
[fsFe,csFe,clFe] = equilibrium(T,cFe,0,TFe1,TFe2,cphsFe1,cphsFe2,...
                                           perTFe,percsFe,perclFe,clap,PhDgFe);
[fsSi,csSi,clSi] = equilibrium(T,cSi,0,TSi1,TSi2,cphsSi1,cphsSi2,...
                                           perTSi,percsSi,perclSi,clap,PhDgSi);

rhosSi = rhosSi0.*(1 - aT.*(T-perTSi) - gCSi.*(csSi-cphsSi1));
rholSi = rholSi0.*(1 - aT.*(T-perTSi) - gCSi.*(clSi-cphsSi1));
rhosFe = rhosFe0.*(1 - aT.*(T-perTFe) - gCFe.*(csFe-cphsFe1));
rholFe = rholFe0.*(1 - aT.*(T-perTFe) - gCFe.*(clFe-cphsFe1));

etasSi = zeros(nzP,nxP) + EtasSi0;
etasFe = zeros(nzP,nxP) + EtasFe0;
Eta = ones(nzP,nxP);

a1      = 1; a2 = 0; a3 = 0; b1 = 1; b2 = 0; b3 = 0;
res = 1e3;
tol = 1e-4;
it = 0;
while res > tol
    fsFei = fsFe; fsSii = fsSi; Pi = Pt;

    if Nx<=10; Pt = mean(mean(Pt(2:end-1,2:end-1))).*ones(size(Pt)); end

    % output crystal fraction and solid and liquid compositions
    [fsFe,csFe,clFe] = equilibrium(T,cFe,Pt,TFe1,TFe2,cphsFe1,cphsFe2,...
                                               perTFe,percsFe,perclFe,clap,PhDgFe);
    [fsSi,csSi,clSi] = equilibrium(T,cSi,Pt,TSi1,TSi2,cphsSi1,cphsSi2,...
                                               perTSi,percsSi,perclSi,clap,PhDgSi);
    fsSi    = min(1,max(0,fsSi)); 
    fFes    = min(1,max(0,fsFe));
    flFe    = 1-fsFe;  
    flSi    = 1-fsSi;
    
    up2date;
    
    rhoRef  = mean(mean(rho(2:end-1,2:end-1)));
    Pt      = rhoRef.*gzP.*ZP + P0;
    T       = Tp .* exp(aT./rhoRef./Cp.*Pt);
    res     = norm(fsFe(:)-fsFei(:),2)./norm(fsFe(:)+TINY,2) ...
            + norm(fsSi(:)-fsSii(:),2)./norm(fsSi(:)+TINY,2) ...
            + norm(Pt(:)  -Pi(:)   ,2)./norm(Pt(:)  +TINY,2);
    it      = it+1;
end

fsFe(xFe==0) = 0;
flFe(xFe==0) = 0;
fsSi(xSi==0) = 0;
flSi(xSi==0) = 0;

Ds = xFe.*fsFe.*dEntrFe + xSi.*fsSi.*dEntrSi;       % bulk entropy change

% set phase [vol] fractions
phisFe = xFe.* fsFe .* rho ./ rhosFe;
philFe = xFe.* flFe .* rho ./ rholFe;
phisSi = xSi.* fsSi .* rho ./ rhosSi;
philSi = xSi.* flSi .* rho ./ rholSi;

% set conserved quantities
FsFe  = rho.*xFe.*fsFe; FlFe  = rho.*xFe.*flFe;             % Solid/liquid iron phase fraction density
FsSi  = rho.*xSi.*fsSi; FlSi  = rho.*xSi.*flSi;             % Solid/liquid silicate phase fraction density
S     = rho.*(Cp.*log(T/T0) - aT./rhoRef.*(Pt-P0) + Ds);    % Bulk entropy density
S0    = rho.*(Cp.*log(T0)   - aT./rhoRef.*(P0)    + Ds); 
slFe  = S./rho - Ds;                                        % phase entropies
slSi  = S./rho - Ds;
ssFe  = slFe + dEntrFe;
ssSi  = slSi + dEntrSi;
XFe   = rho.*xFe; XSi = rho.*xSi;                           % mixture Fe/Si system densities
CFe   = rho.*cFe.*xFe;                                      % mixture Fe component density
CSi   = rho.*cSi.*xSi;                                      % mixture Si component density
XFe   = FsFe + FlFe;
XSi   = FsSi + FlSi;
RHO   = XFe + XSi;                                          % dynamic density

n26Al = 0;
if radheat
% initialise radiogenic isotopes 
Alfrac      = (cSi-cphsSi1)/(cphsSi2-cphsSi1);
nAl0        = nAl_C./mean(xSi(:))./mean(Alfrac(:));    % Bulk silicate Al abundance [kg^{-1}]
n26Al0      = nAl0*Al26_27; % initial numer of 26Al per kg Silicate at CAI formation
tauAl       = t_halfAl/log(2); % lifetime of 26Al
n26Al       = n26Al0*exp(-t_form/tauAl); % initial 26Al content at planetesimal formation

n26Alo      = n26Al;        % initial concentration
dndto       = dndt;         % initial rate of change
end

% reaction transfer rates
GsFe = 0.*fsFe;
GsSi = 0.*fsSi;
GlFe = 0.*flFe;
GlSi = 0.*flSi;


%% initialise previous solution and auxiliary fields
So          = S;
XFeo        = XFe;
XSio        = XSi;
CSio        = CSi;
CFeo        = CFe;
FsFeo       = FsFe;
FsSio       = FsSi;
FlFeo       = FlFe;
FlSio       = FlSi;
dSdto       = dSdt;
dXFedto     = dXFedt;
dXSidto     = dXSidt;
dCFedto     = dCFedt;
dCSidto     = dCSidt;
dFsFedto    = dFsFedt;
dFsSidto    = dFsSidt;
dFlFedto    = dFlFedt;
dFlSidto    = dFlSidt;
rhoo        = rho;
advn_RHOo   = advn_RHO;   
Div_Vo      = Div_V;
dto         = dt;
upd_S       = 0.*dSdt;
upd_XFe     = 0.*dXFedt;
upd_XSi     = 0.*dXSidt;
upd_CFe     = 0.*dCFedt;
upd_CSi     = 0.*dCSidt;
upd_FsFe    = 0.*dFsFedt;
upd_FsSi    = 0.*dFsSidt;
upd_FlFe    = 0.*dFlFedt;
upd_FlSi    = 0.*dFlSidt;
upd_rho     = 0.*dSdt;
dsumMdt     = 0; dsumMdto = dsumMdt;
dsumSdt     = 0; dsumSdto = dsumSdt;
dsumXFedt   = 0; dsumXFedto = dsumXFedt;
dsumXSidt   = 0; dsumXSidto = dsumXSidt;
dsumCFedt   = 0; dsumCFedto = dsumCFedt;
dsumCSidt   = 0; dsumCSidto = dsumCSidt;

%%  additional auxilary fields
hassolSi = flSi<1;
hassolFe = flFe<1;
hasliqSi = fsSi<1;
hasliqFe = fsFe<1;

%% update nonlinear material properties
up2date;
% solve fluid-mechanical equations
switch mode
    case 'spherical'; solve_fluidmech_sp;
    otherwise 
        if Nx==1 
            solve_fluidmech_sp; 
        else 
            solve_fluidmech; 
        end
end
%% check reynold's number
% Velbar = abs(sqrt(((W(1:end-1,2:end-1)+W(2:end,2:end-1))/2).^2 ...
%          + ((U(2:end-1,1:end-1)+U(2:end-1,2:end))/2).^2));
% Re = D*rho(2:end-1,2:end-1).*Vel./Eta(2:end-1,2:end-1);
% 
% figure(45); imagesc(xP(inx),zP(inz),Re); colorbar; axis ij tight


%% initialise recording of model history
switch mode 
    case 'spherical'; history_sp;
    otherwise; history;
end

%% initialise counting variables
RUN.frame   = 0;      % initialise output frame count
time        = 0;      % initialise time count
step        = 0;      % initialise time step count
iter        = 1;      % initialise iteration count
dtlimit     = 'none'; % initialise time limiter

% overwrite fields from file if restarting run
if restart
    if exist(name,'file')
        %path specified in the run script
        fprintf('\n   restart from %s \n\n',name);
        load(name);
        % load(name_h);

        SOL = [W(:);U(:);P(:)];
        RHO = FlFe+FlSi+FsFe+FsSi;
        slFe  = (S - FsFe.*dEntrFe - FsSi.*dEntrSi)./(FlFe+FsFe+FlSi+FsSi);
        slSi  = slFe;
        ssFe  = slFe + dEntrFe;
        ssSi  = slSi + dEntrSi;

        up2date; 

        So          = S;
        XFeo        = XFe;
        XSio        = XSi;
        CSio        = CSi;
        CFeo        = CFe;
        FsFeo       = FsFe;
        FsSio       = FsSi;
        FlFeo       = FlFe;
        FlSio       = FlSi;
        dSdto       = dSdt;
        dXFedto     = dXFedt;
        dXSidto     = dXSidt;
        dCFedto     = dCFedt;
        dCSidto     = dCSidt;
        dFsFedto    = dFsFedt;
        dFsSidto    = dFsSidt;
        dFlFedto    = dFlFedt;
        dFlSidto    = dFlSidt;
        rhoo        = rho;
        advn_RHOo   = advn_RHO;
        Div_Vo      = Div_V;
        dto         = dt;
        dsumMdt     = 0; dsumMdto = dsumMdt;
        dsumSdt     = 0; dsumSdto = dsumSdt;
        dsumXFedt   = 0; dsumXFedto = dsumXFedt;
        dsumXSidt   = 0; dsumXSidto = dsumXSidt;
        dsumCFedt   = 0; dsumCFedto = dsumCFedt;
        dsumCSidt   = 0; dsumCSidto = dsumCSidt;

        %reset history record to current step
        % HST.time = step; if radheat; HST.n26Al = n26Al; end; HST.sumM = HST.sumM(step); HST.sumS = HST.sumS(step);
        % HST.sumXFe  = HST.sumXFe(1:step); HST.sumXSi = HST.sumXSi(1:step); HST.sumCFe = HST.sumCFe(1:step); HST.sumCSi = HST.sumCSi(1:step);
        % HST.dM = HST.dM(1:step); HST.dS = HST.dS(1:step); HST.dXFe = HST.dXFe(1:step); HST.dXSi = HST.dXSi(1:step);
        % HST.dCFe = HST.dCFe(1:step); HST.dCSi = HST.dCSi(1:step); HST.EM = HST.EM(1:step); HST.ES = HST.ES(1:step);
        % HST.EXFe = HST.EXFe(1:step); HST.EXSi = HST.EXSi(1:step); HST.ECFe = HST.ECFe(1:step); HST.ECSi = HST.ECSi(1:step); 

        % HST.time = step; if radheat; HST.n26Al = n26Al; end; HST.sumM = HST.sumM(step); HST.sumS = HST.sumS(step);
        % HST.sumXFe  = HST.sumXFe(1:step); HST.sumXSi = HST.sumXSi(1:step); HST.sumCFe = HST.sumCFe(1:step); HST.sumCSi = HST.sumCSi(1:step);
        % HST.dM = HST.dM(1:step); HST.dS = HST.dS(1:step); HST.dXFe = HST.dXFe(1:step); HST.dXSi = HST.dXSi(1:step);
        % HST.dCFe = HST.dCFe(1:step); HST.dCSi = HST.dCSi(1:step); HST.EM = HST.EM(1:step); HST.ES = HST.ES(1:step);
        % HST.EXFe = HST.EXFe(1:step); HST.EXSi = HST.EXSi(1:step); HST.ECFe = HST.ECFe(1:step); HST.ECSi = HST.ECSi(1:step);

        if radheat
           % n26Al = HST.n26Al(end); % n26Al0*exp(-(t_form+time)/tauAl);
           n26Alo      = n26Al;
        end
      
        % output;

        time    = time+dt;
        step = 1;
        % step    = step+1;

    else % continuation file does not exist, start from scratch
        fprintf('\n   !!! restart file does not exist !!! \n   => starting run from scratch %s \n\n',RunID);
        solve_fluidmech;
        up2date;
        switch mode      
            case 'spherical'; history_sp;     
            otherwise; history; end
        output;
    end
else
    % complete, plot, and save initial condition
    switch mode
        case 'spherical'; solve_fluidmech_sp;
        otherwise; solve_fluidmech;
    end
    up2date;
    switch mode
        case 'spherical'; history_sp;
        otherwise; history;
    end
    output;
    step = step+1;
end