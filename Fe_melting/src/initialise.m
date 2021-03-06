% suspension flow: initialise model run

% print initialisation header
fprintf(1,'  ---  initialise model run \n\n');


%% setup numerical grid
% set dimensions of staggered/ghosted 2D grid
NUM.Nx          =  NUM.N*NUM.L/NUM.D;                   % number of grid points in x-direction
NUM.Nz          =  NUM.N;                               % number of grid points in z-direction
NUM.nxC         =  NUM.Nx+1;                            % number of corner nodes in x-direction
NUM.nzC         =  NUM.Nz+1;                            % number of corner nodes in z-direction
NUM.nxP         =  NUM.Nx+2;                            % number of centre nodes in x-direction
NUM.nzP         =  NUM.Nz+2;                            % number of centre nodes in z-direction
NUM.nxW         =  NUM.Nx+2;                            % number of z-face nodes in x-direction
NUM.nzW         =  NUM.Nz+1;                            % number of z-face nodes in z-direction
NUM.nxU         =  NUM.Nx+1;                            % number of x-face nodes in x-direction
NUM.nzU         =  NUM.Nz+2;                            % number of x-face nodes in z-direction
NUM.NC          =  NUM.nxC*NUM.nzC;                     % total number of corner nodes
NUM.NP          =  NUM.nxP*NUM.nzP;                     % total number of corner nodes
NUM.NW          =  NUM.nxW*NUM.nzW;                     % total number of z-face nodes
NUM.NU          =  NUM.nxU*NUM.nzU;                     % total number of x-face nodes
NUM.NDOF        =  NUM.NP+NUM.NW+NUM.NU;                % total number of all degrees of freedum

% set coordinate vectors
NUM.xC          =  0:NUM.h:NUM.L;                       % Horizontal coordinates of corner nodes [m]
NUM.zC          =  0:NUM.h:NUM.D;                       % Vertical   coordinates of corner nodes [m]
NUM.xP          =  -NUM.h/2:NUM.h:NUM.L+NUM.h/2;        % Horizontal coordinates of centre nodes [m]
NUM.zP          =  -NUM.h/2:NUM.h:NUM.D+NUM.h/2;        % Vertical   coordinates of centre nodes [m]
NUM.xW          =  NUM.xP;                              % Horizontal coordinates of z-face nodes [m]
NUM.zW          =  NUM.zC;                              % Vertical   coordinates of z-face nodes [m]
NUM.xU          =  NUM.xC;                              % Horizontal coordinates of x-face nodes [m]
NUM.zU          =  NUM.zP;                              % Vertical   coordinates of x-face nodes [m]

% set 2D coordinate grids
[NUM.XC,NUM.ZC] = meshgrid(NUM.xC,NUM.zC);              % corner nodes grid
[NUM.XP,NUM.ZP] = meshgrid(NUM.xP,NUM.zP);              % centre nodes grid
[NUM.XW,NUM.ZW] = meshgrid(NUM.xW,NUM.zW);              % z-face nodes grid
[NUM.XU,NUM.ZU] = meshgrid(NUM.xU,NUM.zU);              % x-face nodes grid


%% setup mapping arrays
NUM.Map         =  reshape(1:NUM.NP,NUM.nzP,NUM.nxP);
NUM.MapW        =  reshape(1:NUM.NW,NUM.nzW,NUM.nxW);
NUM.MapU        =  reshape(1:NUM.NU,NUM.nzU,NUM.nxU) + NUM.NW;
NUM.MapP        =  reshape(1:NUM.NP,NUM.nzP,NUM.nxP) + NUM.NW + NUM.NU;

% get smoothed initialisation field
rng(15);
rp = randn(NUM.nzP,NUM.nxP);
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
PHY.gxP = zeros(NUM.nzP,NUM.nxP) + PHY.gx;      PHY.gzP = zeros(NUM.nzP,NUM.nxP) + PHY.gz;
PHY.gx  = zeros(NUM.nzU,NUM.nxU) + PHY.gx;      PHY.gz  = zeros(NUM.nzW,NUM.nxW) + PHY.gz;

      


%% setup velocity-pressure solution arrays
SOL.W           = zeros(NUM.nzW,NUM.nxW);               % z-velocity on z-face nodes
SOL.U           = zeros(NUM.nzU,NUM.nxU);               % x-velocity on x-face nodes
SOL.P           = zeros(NUM.nzP,NUM.nxP);               % pressure on centre nodes
S               = [SOL.W(:);SOL.U(:);SOL.P(:)];         % full solution vector

% projected velocities on centre nodes
SOL.UP          = zeros(NUM.nzP,NUM.nxP);
SOL.WP          = zeros(NUM.nzP,NUM.nxP);
SOL.UP(:,2:end-1) = (SOL.U(:,1:end-1)+SOL.U(:,2:end))./2;
SOL.WP(2:end-1,:) = (SOL.W(1:end-1,:)+SOL.W(2:end,:))./2;

%% setup deformation property arrays
DEF.ups = zeros(NUM.nzP,NUM.nxP);               % velocity divergence on centre nodes
DEF.exx = zeros(NUM.nzP,NUM.nxP);               % x-normal strain rate on centre nodes
DEF.ezz = zeros(NUM.nzP,NUM.nxP);               % z-normal strain rate on centre nodes
DEF.exz = zeros(NUM.nzC,NUM.nxC);               % xz-shear strain rate on corner nodes
DEF.eII = zeros(NUM.nzP,NUM.nxP);               % strain rate magnitude on centre nodes
DEF.txx = zeros(NUM.nzP,NUM.nxP);               % x-normal stress on centre nodes
DEF.tzz = zeros(NUM.nzP,NUM.nxP);               % z-normal stress on centre nodes
DEF.txz = zeros(NUM.nzC,NUM.nxC);               % xz-shear stress on corner nodes
DEF.tII = zeros(NUM.nzP,NUM.nxP);               % stress magnitude on centre nodes


%% setup heating rates
dHdt     = zeros(NUM.Nz+2,NUM.Nx+2);            % enthalpy rate of change
dXdt     = zeros(NUM.Nz+2,NUM.Nx+2);            % Iron system rate of change 
dCSidt   = zeros(NUM.Nz+2,NUM.Nx+2);            % Silicate component density rate of change
dCFedt   = zeros(NUM.Nz+2,NUM.Nx+2);            % Iron component density rate of change
dfFedt   = zeros(NUM.Nz+2,NUM.Nx+2);            % Iron melt fraction rate of change
dfSidt   = zeros(NUM.Nz+2,NUM.Nx+2);            % Silicate melt fraction rate of change
diff_T   = zeros(NUM.Nz+2,NUM.Nx+2);            % Temperature diffusion rate
diff_CSi = zeros(NUM.Nz+2,NUM.Nx+2);            % Silicate component diffusion rate
diff_CFe = zeros(NUM.Nz+2,NUM.Nx+2);            % Iron component diffusion rate
Div_V    = zeros(NUM.Nz+2,NUM.Nx+2);            % Stokes velocity divergence
Div_rhoV = zeros(NUM.Nz+2,NUM.Nx+2);            % Mixture mass flux divergence

%% initialise previous solution and auxiliary fields
rhoo      = ones(NUM.Nz+2,NUM.Nx+2);
Ho        = ones(NUM.Nz+2,NUM.Nx+2);
dHdto     = dHdt;
Div_rhoVo = Div_rhoV;

% initialise counting variables
RUN.frame = 0;      % initialise output frame count
NUM.time  = 0;      % initialise time count
NUM.step  = 0;      % initialise time step count
iter      = 0;      % initialise iteration count

%% set initial condition on thermochemical fields
% temperature
pert = -NUM.h/2.*cos(NUM.XP*2*pi/NUM.D);
% temporary, set theral distribution parameters. Switch to run scripts
% later
zlay     =  0.03;                 % layer thickness (relative to domain depth D)
wlay_T   =  0.05;                % thickness of smooth layer boundary (relative to domain depth D)
% wlay_c   =  2*NUM.h/NUM.D;       % thickness of smooth layer boundary (relative to domain depth D)


switch SOL.Ttype
    case 'constant'     % constant temperature
%         SOL.T      = zeros(NUM.nzP,NUM.nxP) + SOL.T1 + (SOL.T0-SOL.T1).*exp(-NUM.ZP./(10*NUM.h));
        SOL.T      = SOL.T0 + (SOL.T1-SOL.T0) .* (1+erf((NUM.ZP/NUM.D-zlay)/wlay_T))/2; % + dT.*rp;
        SOL.T(1,:) = SOL.T0;
    case 'linear'       % linear temperature gradient with depth
        SOL.T      = SOL.T0 + abs(NUM.ZP+pert)./NUM.D.*(SOL.T1-SOL.T0);
    case 'gaussian'     % constant temperature with gaussian plume
        SOL.T      = zeros(NUM.nzP,NUM.nxP) + SOL.T0;
        SOL.T      = SOL.T + SOL.dT.*exp(-(NUM.XP-SOL.xT).^2./SOL.rT.^2 - (NUM.ZP-SOL.zT).^2./SOL.rT.^2 );
    case 'hot bottom'
        SOL.T      = zeros(NUM.nzP,NUM.nxP) + SOL.T0;
        SOL.T(end-10:end,:) = SOL.T0+100;
%         SOL.T = SOL.T + 1./(1+exp(-(NUM.ZP-SOL.zT+pert)./(NUM.D/50))) .* (SOL.T1-SOL.T0);
end

% set initial component weight fraction [kg/kg]
CHM.xFe = CHM.xFe0 + dxFe.*rp; % Fe system
CHM.xSi = 1 - CHM.xFe;         % Si system
CHM.cFe = CHM.cFe0.*(CHM.xFe>0) + dcFe.*rp; % Fe component
CHM.cSi = CHM.cSi0.*(CHM.xSi>0) + dcSi.*rp; % Si component


% estimate mixture density from Fe/Si system fractions
MAT.rho = CHM.xFe.*PHY.rhoFes + CHM.xSi.*PHY.rhoSis;

res = 1e3;
tol = 1e-16;
it = 0;
while res > tol
    
    rhoRef     = mean(mean(MAT.rho(2:end-1,2:end-1)));
    SOL.Pt     = rhoRef.*PHY.gzP.*NUM.ZP + SOL.P;
    if NUM.Nx<=10; SOL.Pt = mean(mean(SOL.Pt(2:end-1,2:end-1))).*ones(size(SOL.Pt)); end

    % output crystal fraction and fertile solid and liquid concentrations
    [CHM.fFes,CHM.csFe,CHM.clFe] = equilibrium_single(SOL.T,CHM.cFe,SOL.Pt,CHM.TFe1,CHM.TFe2,CHM.cphsFe1,CHM.cphsFe2,...
                                               CHM.perTFe,CHM.perCsFe,CHM.perClFe,CHM.clap,CHM.PhDgFe,TINY);
    [CHM.fSis,CHM.csSi,CHM.clSi] = equilibrium(SOL.T,CHM.cSi,SOL.Pt,CHM.TSi1,CHM.TSi2,CHM.cphsSi1,CHM.cphsSi2,...
                                               CHM.perTSi,CHM.perCsSi,CHM.perClSi,CHM.clap,CHM.PhDgSi,TINY);
    CHM.fFel = 1-CHM.fFes;  CHM.fSil = 1-CHM.fSis;

    up2date;
    
    res        = rhoRef - mean(mean(MAT.rho(2:end-1,2:end-1)));
    it = it+1;
    disp([num2str(it), 'iter'])

end

MAT.phiFes = CHM.xFe.* CHM.fFes .* MAT.rho ./ MAT.rhoFes;
MAT.phiFel = CHM.xFe.* CHM.fFel .* MAT.rho ./ MAT.rhoFel;
MAT.phiSis = CHM.xSi.* CHM.fSis .* MAT.rho ./ MAT.rhoSis;
MAT.phiSil = CHM.xSi.* CHM.fSil .* MAT.rho ./ MAT.rhoSil;

MAT.Ds    = (CHM.xFe.*(CHM.fFes.*0 + CHM.fFel.*CHM.dEntrFe)...        % mixture latent heat capacity density
          +  CHM.xSi.*(CHM.fSis.*0 + CHM.fSil.*CHM.dEntrSi));
    
SOL.H   = SOL.T.*MAT.rho.*(MAT.Ds + PHY.Cp);                                  % mixture enthalpy density
CHM.XFe = MAT.rho.*CHM.xFe; CHM.XSi = MAT.rho.*CHM.xSi;                    % mixture Fe/Si system densities
CHM.CFe = MAT.rho.*CHM.cFe.*CHM.xFe;                                       % mixture Fe component density
CHM.CSi = MAT.rho.*CHM.cSi.*CHM.xSi;                                       % mixture Si component density
CHM.cSi = CHM.CSi./(CHM.xSi+TINY)./MAT.rho;
CHM.cFe = CHM.CFe./(CHM.xFe+TINY)./MAT.rho;

CHM.GFe = 0.*CHM.fFel;
CHM.GSi = 0.*CHM.fSil;


%% initialise radioactive decay (make sure to cleanup)
nAl     = 2.0532e23;% initial nAl per kg
rAl0    = 5.25e-5;  % initial ratio of 26Al/27Al
eAl     = 5e-13;    % decay energy
NAl     = zeros(NUM.nzP,NUM.nxP); 
NAl     = NAl + nAl*rAl0 .* mean(MAT.rho(:)); % initial NAl per m^3
dNdt = zeros(NUM.nzP,NUM.nxP);

MAT.Hr = zeros(NUM.nzP,NUM.nxP) + PHY.Hr0;
HIST.Hr = 0;
if ~RUN.rad; MAT.Hr + PHY.Hr0;  end

%% initialise previous solution and auxiliary fields
rhoo      = MAT.rho;
Ho        = SOL.H;
dHdto     = dHdt;
Div_rhoVo = Div_rhoV;

% initialise counting variables
RUN.frame = 0;      % initialise output frame count
NUM.time  = 0;      % initialise time count
NUM.step  = 0;      % initialise time step count
iter      = 0;      % initialise iteration count


%% update nonlinear material properties
up2date;

%% initialise recording of model history
history;