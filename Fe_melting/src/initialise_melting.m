% suspension flow: initialise model run

% print initialisation header
fprintf(1,'  ---  initialise model run \n\n');


%% setup numerical grid
% set dimensions of staggered/ghosted 2D grid
NUM.nxC         =  NUM.N+1;                             % number of corner nodes in x-direction
NUM.nzC         =  NUM.N+1;                             % number of corner nodes in z-direction
NUM.nxP         =  NUM.N+2;                             % number of centre nodes in x-direction
NUM.nzP         =  NUM.N+2;                             % number of centre nodes in z-direction
NUM.nxW         =  NUM.N+2;                             % number of z-face nodes in x-direction
NUM.nzW         =  NUM.N+1;                             % number of z-face nodes in z-direction
NUM.nxU         =  NUM.N+1;                             % number of x-face nodes in x-direction
NUM.nzU         =  NUM.N+2;                             % number of x-face nodes in z-direction
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
% materials
MAT.Eta         = zeros(NUM.nzP,NUM.nxP) + PHY.Eta0;   	% viscosity
MAT.Eta0        = MAT.Eta;                              % reference viscosity
MAT.aTFe        = zeros(NUM.nzP,NUM.nxP) + PHY.aT0;    	% iron thermal expansivity
MAT.aTSi        = zeros(NUM.nzP,NUM.nxP) + PHY.aT0;    	% silicate thermal expansivity
MAT.kC          = zeros(NUM.nzP,NUM.nxP) + PHY.kC;      % chemical diffusivity
% MAT.kTSi        = zeros(NUM.nzP,NUM.nxP) + PHY.kTSi;    % silicate thermal conductivity
% MAT.kTFe        = zeros(NUM.nzP,NUM.nxP) + PHY.kTFe;    % iron thermal conductivity
% MAT.Cp          = zeros(NUM.nzP,NUM.nxP) + PHY.Cp0;   	% heat capacity
% MAT.Hr	 = zeros(NUM.nzP,NUM.nxP) + PHY.Hr0;	% radiogenic heating
MAT.Hr          = PHY.Hr0 + PHY.Hr0*0.01 .*rp;

%% setup velocity-pressure solution arrays
SOL.W           = zeros(NUM.nzW,NUM.nxW);               % z-velocity on z-face nodes
SOL.U           = zeros(NUM.nzU,NUM.nxU);               % x-velocity on x-face nodes
SOL.P           = zeros(NUM.nzP,NUM.nxP);               % pressure on centre nodes
SOL.Stokes      = zeros(NUM.nzP,NUM.nxP);               % Stokes settling velocity

% projected velocities on centre nodes
SOL.UP          = zeros(NUM.nzP,NUM.nxP);
SOL.WP          = zeros(NUM.nzP,NUM.nxP);
SOL.UP(:,2:end-1) = SOL.U(:,1:end-1)+SOL.U(:,2:end)./2;
SOL.WP(2:end-1,:) = SOL.W(1:end-1,:)+SOL.W(2:end,:)./2;
S               = [SOL.W(:);SOL.U(:);SOL.P(:)];             % Solution matrix

%% set thermochemical fields
% temperature
pert = -NUM.h/2.*cos(NUM.XP*2*pi/NUM.D);
switch SOL.Ttype
    case 'constant'     % constant temperature
        SOL.T  = zeros(NUM.nzP,NUM.nxP) + SOL.T0;
    case 'linear'       % linear temperature gradient with depth
        SOL.T  = SOL.T0 + abs(NUM.ZP+pert)./NUM.D.*(SOL.T1-SOL.T0);
    case 'gaussian'     % constant temperature with gaussian plume
        SOL.T = zeros(NUM.nzP,NUM.nxP) + SOL.T0;
        SOL.T = SOL.T + SOL.dT.*exp(-(NUM.XP-SOL.xT).^2./SOL.rT.^2 - (NUM.ZP-SOL.zT).^2./SOL.rT.^2 );
    case 'hot bottom'
        SOL.T = zeros(NUM.nzP,NUM.nxP) + SOL.T0;
        SOL.T(end-10:end,:) = SOL.T0+100;
%         SOL.T = SOL.T + 1./(1+exp(-(NUM.ZP-SOL.zT+pert)./(NUM.D/50))) .* (SOL.T1-SOL.T0);
end
% convert from potential to natural temperature (not important yet, will
% consider reintroducing later)
% SOL.T            = SOL.T.*exp(PHY.aT0*(PHY.gzP.*(NUM.ZP+pert) + PHY.gxP.*NUM.XP)./PHY.Cp0);
% SOL.T([1 end],:) = SOL.T([2 end-1],:);
% SOL.T(:,[1 end]) = SOL.T(:,[2 end-1]);

% set initial component weight fraction [kg/kg]
cFe                 = cFe0 + dc.*rp; % Fe component
cSi                 = cSi0 + dc.*rp; % Si component
xFe                 = xFe0 + dc.*rp; % Fe system
xSi                 = 1 - xFe;
% estimate pressure from estimated wt%
SOL.Pt = (xFe.*PHY.rhoFes + xSi.*PHY.rhoSis).*NUM.ZP.*PHY.gzP;

% output crystal fraction and fertile solid and liquid concentrations
[fFes,csFe,clFe]     = equilibrium(SOL.T,cFe,SOL.Pt,TFe1,TFe2,cphsFe1,cphsFe2,...
                                  perTSi,perCsFe,perClFe,clap,PhDgFe,TINY);
[fSis,csSi,clSi]     = equilibrium(SOL.T,cSi,SOL.Pt,TSi1,TSi2,cphsSi1,cphsSi2,...
                                  perTSi,perCsSi,perClSi,clap,PhDgSi,TINY);
fFel = 1-fFes; fSil = 1-fSis;

MAT.rhoSis = PHY.rhoSis.*(1 - MAT.aTSi.*(SOL.T-SOL.T0) - PHY.gammaSi.*cSi);
MAT.rhoSil = PHY.rhoSil.*(1 - MAT.aTSi.*(SOL.T-SOL.T0) - PHY.gammaSi.*cSi);
MAT.rhoFes = PHY.rhoFes.*(1 - MAT.aTFe.*(SOL.T-SOL.T0) - PHY.gammaFe.*cFe);
MAT.rhoFel = PHY.rhoFel.*(1 - MAT.aTFe.*(SOL.T-SOL.T0) - PHY.gammaFe.*cFe);

MAT.rhot   = 1./ (xFe.*fFes./MAT.rhoFes + xFe.*fFel./MAT.rhoFel + xSi.*fSis./MAT.rhoSis + xSi.*fSil./MAT.rhoSil);

phiFes     = xFe.* fFes .* MAT.rhot ./ MAT.rhoFes;
phiFel     = xFe.* fFel .* MAT.rhot ./ MAT.rhoFel;
phiSis     = xSi.* fSis .* MAT.rhot ./ MAT.rhoSis;
phiSil     = xSi.* fSil .* MAT.rhot ./ MAT.rhoSil; 

XFe = MAT.rhot.*xFe; XSi = MAT.rhot.*xSi;
 SOL.Pt = (XFe.*(fFes.*PHY.rhoFes + fFel.*PHY.rhoFel)...
        +  XSi.*(fSis.*PHY.rhoSis + fSil.*PHY.rhoSil)).*NUM.ZP.*PHY.gzP;


% % loop for P-dependent equilibrium to be added later on
% % start loop to solve P, T and phi nonlinearities
% % initial boussineq approximation
% phiFes = xFe.*fFes; phiFel = xFe.*fFel; phiSis = xSi.*fSis; phiSil = xSi.*fSil;
% MAT.rhot = MAT.rhoSil.*phiSil + MAT.rhoSis.*phiSis + MAT.rhoFel.*phiFel +MAT.rhoFes.*phiFes;
% for i = 1:10
% [fFes,csFe,clFe]     = equilibrium(SOL.T,cFe,SOL.Pt,TFe1,TFe2,cphsFe1,cphsFe2,...
%                                   perTSi,perCsFe,perClFe,clap,PhDgFe,TINY);
% [fSis,csSi,clSi]     = equilibrium(SOL.T,cSi,SOL.Pt,TSi1,TSi2,cphsSi1,cphsSi2,...
%                                   perTSi,perCsSi,perClSi,clap,PhDgSi,TINY);
% fFel = 1-fFes; fSil = 1-fSis;
% 
% phiFesi = phiFes; phiFeli = phiFel; phiSisi = phiSis; phiSili = phiSil;
% phiFes = max(TINY,min(1-TINY, xFe.*fFes.*MAT.rhot./MAT.rhoFes ))./(phiFes+phiFel+phiSis+phiSil);
% phiFel = max(TINY,min(1-TINY, xFe.*fFel.*MAT.rhot./MAT.rhoFel ))./(phiFes+phiFel+phiSis+phiSil);
% phiSis = max(TINY,min(1-TINY, xSi.*fSis.*MAT.rhot./MAT.rhoSis ))./(phiFes+phiFel+phiSis+phiSil);
% phiSil = max(TINY,min(1-TINY, xSi.*fSil.*MAT.rhot./MAT.rhoSil ))./(phiFes+phiFel+phiSis+phiSil);
% 
% MAT.rhot                = phiFes.*MAT.rhoFes...
%                         + phiFel.*MAT.rhoFel...
%                         + phiSis.*MAT.rhoSis...
%                         + phiSil.*MAT.rhoSil;
%                     
% XFe = MAT.rhot.*xFe; XSi = MAT.rhot.*xSi;                    
%                     
% SOL.Pt = (XFe.*(fFes.*PHY.rhoFes + fFel.*PHY.rhoFel)...
%        +  XSi.*(fSis.*PHY.rhoSis + fSil.*PHY.rhoSil)).*NUM.ZP.*PHY.gzP;
%    
% 
% end

rhoref  = mean(mean(MAT.rhot(2:end-1,2:end-1)));
rhoCpt  = (XFe.*(fFes.*PHY.CpFes + fFel.*PHY.CpFel)...
        +  XSi.*(fSis.*PHY.CpSis + fSil.*PHY.CpSil));
SOL.H   =  SOL.T.*(XFe.*(fFel.*dEntrSi) + XSi.*(fSil.*dEntrSi) +rhoCpt);
CFe     =  cFe.*XFe.*MAT.rhot;   CSi     =  cSi.*XSi.*MAT.rhot;   % component densities


%% setup deformation property arrays
DEF.ups         = zeros(NUM.nzP,NUM.nxP);               % velocity divergence on centre nodes
DEF.exx         = zeros(NUM.nzP,NUM.nxP);               % x-normal strain rate on centre nodes
DEF.ezz         = zeros(NUM.nzP,NUM.nxP);               % z-normal strain rate on centre nodes
DEF.exz         = zeros(NUM.nzC,NUM.nxC);               % xz-shear strain rate on corner nodes
DEF.eII         = zeros(NUM.nzP,NUM.nxP);               % strain rate magnitude on centre nodes
DEF.txx         = zeros(NUM.nzP,NUM.nxP);               % x-normal stress on centre nodes
DEF.tzz         = zeros(NUM.nzP,NUM.nxP);               % z-normal stress on centre nodes
DEF.txz         = zeros(NUM.nzC,NUM.nxC);               % xz-shear stress on corner nodes
DEF.tII         = zeros(NUM.nzP,NUM.nxP);               % stress magnitude on centre nodes


%% setup heating rates
% SOL.dTdt        = zeros(NUM.N  ,NUM.N );              % temperature rate of change
dHdt        = zeros(NUM.N+2,NUM.N+2);                % enthalpy rate of change
dXdt        = zeros(NUM.N+2,NUM.N+2); 
dCSidt      = zeros(NUM.N+2,NUM.N+2);                % Silicate component density rate of change
dCFedt      = zeros(NUM.N+2,NUM.N+2);                % Iron component density rate of change
diff_T      = zeros(NUM.N+2,NUM.N+2);
diff_CSi    = zeros(NUM.N+2,NUM.N+2);
diff_CFe    = zeros(NUM.N+2,NUM.N+2);
% SOL.Hs          = zeros(NUM.nzP,NUM.nxP);          	% shear heating rate
% SOL.Ha          = zeros(NUM.nzP,NUM.nxP);          	% adiabatic heating rate
Div_V         = zeros(NUM.N+2  ,NUM.N+2 );              % level set free surface change

% store previous solutions
    MAT.rhoo = MAT.rhot;
    Ho                  = SOL.H;     % store previous temperature solution
    dHdto               = dHdt;  % store previous rate of change
    Div_rhoVo = 0;
%     phio                = SOL.phi;
%% update nonlinear material properties
RUN.frame = 0;     % initialise output frame count
NUM.time  = 0;      % initialise time count
NUM.step  = 0;      % initialise time step count
up2date;


%% output initial condition
% output;

