% suspension flow: initialise model run

% print initialisation header
fprintf(1,'  ---  initialise model run \n\n');


%% setup numerical grid
% set dimensions of staggered/ghosted 2D grid
Nx          =  N*L/D;                           % number of grid points in x-direction
Nz          =  N;                               % number of grid points in z-direction
nxC         =  Nx+1;                            % number of corner nodes in x-direction
nzC         =  Nz+1;                            % number of corner nodes in z-direction
nxP         =  Nx+2;                            % number of centre nodes in x-direction
nzP         =  Nz+2;                            % number of centre nodes in z-direction
nxW         =  Nx+2;                            % number of z-face nodes in x-direction
nzW         =  Nz+1;                            % number of z-face nodes in z-direction
nxU         =  Nx+1;                            % number of x-face nodes in x-direction
nzU         =  Nz+2;                            % number of x-face nodes in z-direction
NC          =  nxC*nzC;                         % total number of corner nodes
NP          =  nxP*nzP;                         % total number of corner nodes
NW          =  nxW*nzW;                         % total number of z-face nodes
NU          =  nxU*nzU;                         % total number of x-face nodes
NDOF        =  NP+NW+NU;                        % total number of all degrees of freedum

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
Map         =  reshape(1:NP,nzP,nxP);
MapW        =  reshape(1:NW,nzW,nxW);
MapU        =  reshape(1:NU,nzU,nxU) + NW;
MapP        =  reshape(1:NP,nzP,nxP);

inz = 2:length(zP)-1;
inx = 2:length(xP)-1;

% get smoothed initialisation field
rng(67);
rp = randn(nzP,nxP);
for i = 1:round(smth)
    rp(2:end-1,2:end-1) = rp(2:end-1,2:end-1) + diff(rp(:,2:end-1),2,1)./8 ...
                                              + diff(rp(2:end-1,:),2,2)./8;
    rp([1 end],:)       = rp([2 end-1],:);
    rp(:,[1 end])       = rp(:,[2 end-1]);
end
rp              = rp./max(abs(rp(:)));
rp              = rp - mean(mean(rp(2:end-1,2:end-1)));


%% setup material property arrays
% gravity
gxP = zeros(nzP,nxP) + gx0;      gzP = zeros(nzP,nxP) + gz0;
gx  = zeros(nzU,nxU) + gx0;      gz  = zeros(nzW,nxW) + gz0;


%% setup velocity-pressure solution arrays
W           = zeros(nzW,nxW);               % z-velocity on z-face nodes
U           = zeros(nzU,nxU);               % x-velocity on x-face nodes
P           = zeros(nzP,nxP);               % pressure on centre nodes
SOL         = [W(:);U(:);P(:)];             % full solution vector
Vel         = zeros(nzP,nxP);               % velocity magnitude

% projected velocities on centre nodes
UP              = zeros(nzP,nxP);
WP              = zeros(nzP,nxP);
UP(:,2:end-1)   = (U(:,1:end-1)+U(:,2:end))./2;
WP(2:end-1,:)   = (W(1:end-1,:)+W(2:end,:))./2;

%% setup deformation property arrays
SOLups = zeros(nzP,nxP);               % velocity divergence on centre nodes
SOLexx = zeros(nzP,nxP);               % x-normal strain rate on centre nodes
SOLezz = zeros(nzP,nxP);               % z-normal strain rate on centre nodes
SOLexz = zeros(nzC,nxC);               % xz-shear strain rate on corner nodes
SOLeII = zeros(nzP,nxP);               % strain rate magnitude on centre nodes
SOLtxx = zeros(nzP,nxP);               % x-normal stress on centre nodes
SOLtzz = zeros(nzP,nxP);               % z-normal stress on centre nodes
SOLtxz = zeros(nzC,nxC);               % xz-shear stress on corner nodes
SOLtII = zeros(nzP,nxP);               % stress magnitude on centre nodes


%% setup heating rates
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
Div_rhoV    = zeros(Nz,Nx);             % Mixture mass flux divergence
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
xFe = xFe0 + dxFe.*rp; % Fe system
xSi = 1 - xFe;         % Si system
cFe = zeros(size(xFe)) + cFe0 + dcFe.*rp; % Fe component
cSi = zeros(size(xSi)) + cSi0 + dcSi.*rp - cSimin; % Si component
Pt = 0;

%% initialise radioactive decay (make sure to cleanup)
% nAl     = 2.0532e23;% initial nAl per kg
% rAl0    = 5.25e-5;  % initial ratio of 26Al/27Al
% eAl     = 5e-13;    % decay energy
% NAl     = zeros(nzP,nxP); 
% NAl     = NAl + nAl*rAl0 .* mean(rho(:)); % initial NAl per m^3
% dNdt = zeros(nzP,nxP);

Hr  = Hr0.*ones(size(T(2:end-1,2:end-1)));

%% initialise loop
[fsFe,csFe,clFe] = equilibrium(T,cFe,0,TFe1,TFe2,cphsFe1,cphsFe2,...
                                           perTFe,percsFe,perclFe,clap,PhDgFe);
[fsSi,csSi,clSi] = equilibrium(T,cSi,0,TSi1,TSi2,cphsSi1,cphsSi2,...
                                           perTSi,percsSi,perclSi,clap,PhDgSi);

rhosSi = rhosSi0.*(1 - aT.*(T-perTSi) - gCSi.*(csSi-cphsSi1));
rholSi = rholSi0.*(1 - aT.*(T-perTSi) - gCSi.*(clSi-cphsSi1));
rhosFe = rhosFe0.*(1 - aT.*(T-perTFe) - gCFe.*(csFe-cphsFe1));
rholFe = rholFe0.*(1 - aT.*(T-perTFe) - gCFe.*(clFe-cphsFe1));

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
RHO   = XFe + XSi;                                          % dynamic density

% reaction transfer rates
GFes = 0.*fsFe;
GSis = 0.*fsSi;
GFel = 0.*flFe;
GSil = 0.*flSi;


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
Div_rhoVo   = Div_rhoV;
Div_Vo      = Div_V;
dto         = dt;
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
solve_fluidmech;

%% check reynold's number
% Velbar = abs(sqrt(((W(1:end-1,2:end-1)+W(2:end,2:end-1))/2).^2 ...
%          + ((U(2:end-1,1:end-1)+U(2:end-1,2:end))/2).^2));
% Re = D*rho(2:end-1,2:end-1).*Vel./Eta(2:end-1,2:end-1);
% 
% figure(45); imagesc(xP(inx),zP(inz),Re); colorbar; axis ij tight


%% initialise recording of model history
history;

%% initialise counting variables
RUN.frame   = 0;      % initialise output frame count
time        = 0;      % initialise time count
step        = 1;      % initialise time step count
iter        = 1;      % initialise iteration count
dtlimit     = 'none'; % initialise time limiter