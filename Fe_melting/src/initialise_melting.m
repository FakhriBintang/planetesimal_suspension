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
        SOL.T  = zeros(NUM.nzP,NUM.nxP) + SOL.T1;
        SOL.T(1,:) = SOL.T0;
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
CHM.cFe                 = CHM.cFe0 + dc.*rp; % Fe component
CHM.cSi                 = CHM.cSi0 + dc.*rp; % Si component
CHM.xFe                 = CHM.xFe0 + dc.*rp; % Fe system
% CHM.xFe(1:15,:)         = CHM.xFe(1:15,:) +0.05
CHM.xSi                 = 1 - CHM.xFe;

% estimate pressure from estimated wt%
SOL.Pt = (CHM.xFe.*PHY.rhoFes + CHM.xSi.*PHY.rhoSis).*NUM.ZP.*PHY.gzP;

% output crystal fraction and fertile solid and liquid concentrations
[CHM.fFes,CHM.csFe,CHM.clFe]     = equilibrium(SOL.T,CHM.cFe,SOL.Pt,CHM.TFe1,CHM.TFe2,CHM.cphsFe1,CHM.cphsFe2,...
                                  CHM.perTFe,CHM.perCsFe,CHM.perClFe,CHM.clap,CHM.PhDgFe,TINY);
[CHM.fSis,CHM.csSi,CHM.clSi]     = equilibrium(SOL.T,CHM.cSi,SOL.Pt,CHM.TSi1,CHM.TSi2,CHM.cphsSi1,CHM.cphsSi2,...
                                  CHM.perTSi,CHM.perCsSi,CHM.perClSi,CHM.clap,CHM.PhDgSi,TINY);
CHM.fFel = 1-CHM.fFes; CHM.fSil = 1-CHM.fSis;

MAT.rhoSis = PHY.rhoSis.*(1 - MAT.aTSi.*(SOL.T-SOL.T0) - PHY.gammaSi.*(CHM.csSi-CHM.cphsSi1));  
MAT.rhoSil = PHY.rhoSil.*(1 - MAT.aTSi.*(SOL.T-SOL.T0) - PHY.gammaSi.*(CHM.clSi-CHM.cphsSi1));
MAT.rhoFes = PHY.rhoFes.*(1 - MAT.aTFe.*(SOL.T-SOL.T0) - PHY.gammaFe.*(CHM.csFe-CHM.cphsFe1));  
MAT.rhoFel = PHY.rhoFel.*(1 - MAT.aTFe.*(SOL.T-SOL.T0) - PHY.gammaFe.*(CHM.clFe-CHM.cphsFe1));

MAT.rhot   = 1./ (CHM.xFe.*CHM.fFes./MAT.rhoFes + CHM.xFe.*CHM.fFel./MAT.rhoFel + CHM.xSi.*CHM.fSis./MAT.rhoSis + CHM.xSi.*CHM.fSil./MAT.rhoSil);

SOL.phiFes     = CHM.xFe.* CHM.fFes .* MAT.rhot ./ MAT.rhoFes;
SOL.phiFel     = CHM.xFe.* CHM.fFel .* MAT.rhot ./ MAT.rhoFel;
SOL.phiSis     = CHM.xSi.* CHM.fSis .* MAT.rhot ./ MAT.rhoSis;
SOL.phiSil     = CHM.xSi.* CHM.fSil .* MAT.rhot ./ MAT.rhoSil; 

CHM.XFe     = MAT.rhot.*CHM.xFe; CHM.XSi = MAT.rhot.*CHM.xSi;
rhoRef  = mean(mean(MAT.rhot(2:end-1,2:end-1)));
SOL.Pt  = rhoRef.*PHY.gzP.*NUM.ZP + SOL.P;

% fh8 = figure(8); clf;
% TT = linspace(CHM.TSi1,CHM.TSi2,1e3);
% cc = [linspace(CHM.cphsSi2,(CHM.perCsSi+CHM.perClSi)/2,round((CHM.perTSi-CHM.TSi1)./(CHM.TSi2-CHM.TSi1)*1e3)),linspace((CHM.perCsSi+CHM.perClSi)/2,CHM.cphsSi1,round((CHM.perTSi-CHM.TSi2)./(CHM.TSi1-CHM.TSi2)*1e3))];
% [~,CCxSi,CClSi]     = equilibrium(TT,cc,0.*TT,CHM.TSi1,CHM.TSi2,CHM.cphsSi1,CHM.cphsSi2,...
%                                   CHM.perTSi,CHM.perCsSi,CHM.perClSi,CHM.clap,CHM.PhDgSi,TINY);
% 
% TT2 = linspace(CHM.TFe1,CHM.TFe2,1000);
% cc2 = linspace(CHM.cphsFe2,CHM.cphsFe1,length(TT2));
% % cc2 = [linspace(CHM.cphsFe2,(CHM.perCsFe+CHM.perClFe)/2,round((CHM.perTFe-CHM.TFe1)./(CHM.TFe2-CHM.TFe1)*1e3)),linspace((CHM.perCsFe+CHM.perClFe)/2,CHM.cphsFe1,round((CHM.perTFe-CHM.TFe2)./(CHM.TFe1-CHM.TFe2)*1e3))];
% [~,CCxFe,CClFe]     = equilibrium(TT2,cc2,0.*TT2,CHM.TFe1,CHM.TFe2,CHM.cphsFe1,CHM.cphsFe2,...
%                                   CHM.perTFe,CHM.perCsFe,CHM.perClFe,CHM.clap,CHM.PhDgFe,TINY);
% subplot(1,2,1)
% plot(CCxSi,TT,'k-','LineWidth',2); axis tight; hold on; box on;
% plot(CClSi,TT,'k-','LineWidth',2);
% 
% subplot(1,2,2)
% plot(CCxFe,TT2,'r-','LineWidth',2); axis tight; hold on; box on;
% plot(CClFe,TT2,'b-','LineWidth',2);
% % loop for P-dependent equilibrium to be added later on
% % start loop to solve P, T and phi nonlinearities
% % initial boussineq approximation
% SOL.phiFes = CHM.xFe.*CHM.fFes; SOL.phiFel = CHM.xFe.*CHM.fFel; SOL.phiSis = CHM.xSi.*CHM.fSis; SOL.phiSil = CHM.xSi.*CHM.fSil;
% MAT.rhot = MAT.rhoSil.*SOL.phiSil + MAT.rhoSis.*SOL.phiSis + MAT.rhoFel.*SOL.phiFel +MAT.rhoFes.*SOL.phiFes;
% for i = 1:10
% [CHM.fFes,CHM.csFe,CHM.clFe]     = equilibrium(SOL.T,CHM.cFe,SOL.Pt,CHM.TFe1,CHM.TFe2,CHM.cphsFe1,CHM.cphsFe2,...
%                                   CHM.perTSi,CHM.perCsFe,CHM.perClFe,CHM.clap,CHM.PhDgFe,TINY);
% [CHM.fSis,CHM.csSi,CHM.clSi]     = equilibrium(SOL.T,CHM.cSi,SOL.Pt,CHM.TSi1,CHM.TSi2,CHM.cphsSi1,CHM.cphsSi2,...
%                                   CHM.perTSi,CHM.perCsSi,CHM.perClSi,CHM.clap,CHM.PhDgSi,TINY);
% CHM.fFel = 1-CHM.fFes; CHM.fSil = 1-CHM.fSis;
% 
% SOL.phiFesi = SOL.phiFes; SOL.phiFeli = SOL.phiFel; SOL.phiSisi = SOL.phiSis; SOL.phiSili = SOL.phiSil;
% SOL.phiFes = max(TINY,min(1-TINY, CHM.xFe.*CHM.fFes.*MAT.rhot./MAT.rhoFes ))./(SOL.phiFes+SOL.phiFel+SOL.phiSis+SOL.phiSil);
% SOL.phiFel = max(TINY,min(1-TINY, CHM.xFe.*CHM.fFel.*MAT.rhot./MAT.rhoFel ))./(SOL.phiFes+SOL.phiFel+SOL.phiSis+SOL.phiSil);
% SOL.phiSis = max(TINY,min(1-TINY, CHM.xSi.*CHM.fSis.*MAT.rhot./MAT.rhoSis ))./(SOL.phiFes+SOL.phiFel+SOL.phiSis+SOL.phiSil);
% SOL.phiSil = max(TINY,min(1-TINY, CHM.xSi.*CHM.fSil.*MAT.rhot./MAT.rhoSil ))./(SOL.phiFes+SOL.phiFel+SOL.phiSis+SOL.phiSil);
% 
% MAT.rhot                = SOL.phiFes.*MAT.rhoFes...
%                         + SOL.phiFel.*MAT.rhoFel...
%                         + SOL.phiSis.*MAT.rhoSis...
%                         + SOL.phiSil.*MAT.rhoSil;
%                     
% CHM.XFe = MAT.rhot.*CHM.xFe; CHM.XSi = MAT.rhot.*CHM.xSi;                    
%                     
% SOL.Pt = (CHM.XFe.*(CHM.fFes.*PHY.rhoFes + CHM.fFel.*PHY.rhoFel)...
%        +  CHM.XSi.*(CHM.fSis.*PHY.rhoSis + CHM.fSil.*PHY.rhoSil)).*NUM.ZP.*PHY.gzP;
%    
% 
% end

rhoCpt  = (MAT.rhot.*CHM.xFe.*(CHM.fFes.*PHY.CpFes + CHM.fFel.*PHY.CpFel)...
        +  MAT.rhot.*CHM.xSi.*(CHM.fSis.*PHY.CpSis + CHM.fSil.*PHY.CpSil));
SOL.H   =  SOL.T.*(MAT.rhot.*CHM.xFe.*(CHM.fFel.*CHM.dEntrSi) + MAT.rhot.*CHM.xSi.*(CHM.fSil.*CHM.dEntrSi) +rhoCpt);
CHM.CFe     =  MAT.rhot.*CHM.cFe.*CHM.xFe;   
CHM.CSi     =  MAT.rhot.*CHM.cSi.*CHM.xSi;   % component densities


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

Mass0       = sum(sum(MAT.rhot(2:end-1,2:end-1)));
MassErr     = 0;
sumH0       = sum(sum(SOL.H(2:end-1,2:end-1)));

%% output initial condition
% output;

