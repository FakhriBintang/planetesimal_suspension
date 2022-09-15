% print update header
% fprintf(1,'  ---  update materials & deformation \n');
tic;  % start clock on update

%% convert weight to volume fraction
% update T and chemical-dependent density
MAT.rhoSis = PHY.rhoSis.*(1 - PHY.aT.*(SOL.T-CHM.perTSi) - PHY.gCSi.*(CHM.csSi-CHM.cphsSi1));
MAT.rhoSil = PHY.rhoSil.*(1 - PHY.aT.*(SOL.T-CHM.perTSi) - PHY.gCSi.*(CHM.clSi-CHM.cphsSi1));
MAT.rhoFes = PHY.rhoFes.*(1 - PHY.aT.*(SOL.T-CHM.TFe1)   - PHY.gCFe.*(CHM.csFe-CHM.cphsFe1));
MAT.rhoFel = PHY.rhoFel.*(1 - PHY.aT.*(SOL.T-CHM.TFe1)   - PHY.gCFe.*(CHM.clFe-CHM.cphsFe1));

% update mixture density
MAT.rho    = 1./(CHM.xFe.*CHM.fFes./MAT.rhoFes + CHM.xFe.*CHM.fFel./MAT.rhoFel ...
               + CHM.xSi.*CHM.fSis./MAT.rhoSis + CHM.xSi.*CHM.fSil./MAT.rhoSil);
MAT.rho([1 end],:) = MAT.rho([2 end-1],:);  MAT.rho(:,[1 end]) = MAT.rho(:,[2 end-1]);

% update volume fractions
MAT.phiFes = CHM.xFe.* CHM.fFes .* MAT.rho ./ MAT.rhoFes;
MAT.phiFel = CHM.xFe.* CHM.fFel .* MAT.rho ./ MAT.rhoFel;
MAT.phiSis = CHM.xSi.* CHM.fSis .* MAT.rho ./ MAT.rhoSis;
MAT.phiSil = CHM.xSi.* CHM.fSil .* MAT.rho ./ MAT.rhoSil; 


%% update viscosity
% get pure phase densities
etas   = zeros(size(MAT.phiSis)) + PHY.EtaSol0; 
etaSil = PHY.EtaSil0 .* exp(Em./(8.3145.*(SOL.T+273.15))-Em./(8.3145.*((CHM.TFe1+CHM.TFe2)/2+273.15)));
etaFel = zeros(size(MAT.phiFel)) + PHY.EtaFel0;

% get permission weights
kv = permute(cat(3,etas,etaSil,etaFel),[3,1,2]);
Mv = permute(repmat(kv,1,1,1,3),[4,1,2,3])./permute(repmat(kv,1,1,1,3),[1,4,2,3]);
 
ff = max(1e-6,min(1-1e-6,permute(cat(3,MAT.phiSis+MAT.phiFes,MAT.phiSil,MAT.phiFel),[3,1,2])));
FF = permute(repmat(ff,1,1,1,3),[4,1,2,3]);
Sf = (FF./BBP).^(1./CCP);  Sf = Sf./sum(Sf,2);
Xf = sum(AAP.*Sf,2).*FF + (1-sum(AAP.*Sf,2)).*Sf;

% get momentum and volume permissions
thtv = squeeze(prod(Mv.^Xf,2));

% get momentum flux and transfer coefficients
Cv = ((1-ff)./[PHY.dx;.1e-16;PHY.df].^2.*kv.*thtv).^-1;

% compose effective viscosity, segregation coefficients
MAT.Eta  = squeeze(sum(ff.*kv.*thtv,1));                                             % effective magma viscosity
MAT.Eta  = (1./etamax + 1./(MAT.Eta + etamin)).^(-1);                       % limit viscosity range
MAT.Eta([1 end],:) = MAT.Eta([2 end-1],:);  
MAT.Eta(:,[1 end]) = MAT.Eta(:,[2 end-1]);
MAT.EtaC = (MAT.Eta(1:end-1,1:end-1)+MAT.Eta(2:end,1:end-1) ...            % viscosity in cell corners
         +  MAT.Eta(1:end-1,2:end  )+MAT.Eta(2:end,2:end  ))./4;

% Ksgr_x   = max(1e-15,min(1e-6,(MAT.phiSis+MAT.phiFes)./squeeze(Cv(1,:,:))));
% Ksgr_m   = max(1e-15,min(1e-6, MAT.phiSil            ./squeeze(Cv(2,:,:))));
% Ksgr_f   = max(1e-15,min(1e-6, MAT.phiFel            ./squeeze(Cv(3,:,:))));
%test alternative
Ksgr_x = squeeze(Cv(1,:,:)) + 1e-18;  Ksgr_x([1 end],:) = Ksgr_x([2 end-1],:);  Ksgr_x(:,[1 end]) = Ksgr_x(:,[2 end-1]);
Ksgr_m = squeeze(Cv(2,:,:)) + 1e-18;  Ksgr_m([1 end],:) = Ksgr_m([2 end-1],:);  Ksgr_m(:,[1 end]) = Ksgr_m(:,[2 end-1]);
Ksgr_f = squeeze(Cv(3,:,:)) + 1e-18;  Ksgr_f([1 end],:) = Ksgr_f([2 end-1],:);  Ksgr_f(:,[1 end]) = Ksgr_f(:,[2 end-1]);




%% other bulk thermochemical properties
% Update pressure
SOL.Pt      = rhoRef.*MAT.gzP.*NUM.ZP + SOL.P + P0;

% update heat capacities
% MAT.rhoCp = (MAT.rho.*CHM.xFe.*(CHM.fFes.*PHY.CpFes + CHM.fFel.*PHY.CpFel)...  % mixture sensible heat capacity density
%           +  MAT.rho.*CHM.xSi.*(CHM.fSis.*PHY.CpSis + CHM.fSil.*PHY.CpSil));
% MAT.rhoDs = (MAT.rho.*CHM.xFe.*CHM.fFel.*CHM.dEntrFe...                        % mixture latent heat capacity density
%           +  MAT.rho.*CHM.xSi.*CHM.fSil.*CHM.dEntrSi);
MAT.Ds    =  CHM.xFe.*CHM.fFel.*CHM.dEntrFe + CHM.xSi.*CHM.fSil.*CHM.dEntrSi;   % mixture entropy
MAT.ks    =  (CHM.xFe.*PHY.kTFe + CHM.xSi.*PHY.kTSi)./SOL.T;                            % magma thermal conductivity

%% calculate stokes settling velocity
sds = SOL.BCsides;      % side boundary type 
% also runge kutta advection time stepping
kappaseg = 100;

% test alternative
segSis = (((MAT.rhoSis(1:end-1,:)+MAT.rhoSis(2:end,:))./2-(MAT.rho(1:end-1,:)+MAT.rho(2:end,:))./2).*PHY.gz).*2./(1./Ksgr_x(1:end-1,:)+1./Ksgr_x(2:end,:));
segFes = (((MAT.rhoFes(1:end-1,:)+MAT.rhoFes(2:end,:))./2-(MAT.rho(1:end-1,:)+MAT.rho(2:end,:))./2).*PHY.gz).*2./(1./Ksgr_x(1:end-1,:)+1./Ksgr_x(2:end,:));
segSil = (((MAT.rhoSil(1:end-1,:)+MAT.rhoSil(2:end,:))./2-(MAT.rho(1:end-1,:)+MAT.rho(2:end,:))./2).*PHY.gz).*2./(1./Ksgr_m(1:end-1,:)+1./Ksgr_m(2:end,:));
segFel = (((MAT.rhoFel(1:end-1,:)+MAT.rhoFel(2:end,:))./2-(MAT.rho(1:end-1,:)+MAT.rho(2:end,:))./2).*PHY.gz).*2./(1./Ksgr_f(1:end-1,:)+1./Ksgr_f(2:end,:));


% zero segregation at the boundaries for accumulation
for nsmooth = 1:nvsmooth
    segSis(2:end-1,2:end-1)   = segSis(2:end-1,2:end-1) + diff(segSis(:,2:end-1),2,1)./8 + diff(segSis(2:end-1,:),2,2)./8;
    segSis([1 end],:) = 0;
    segSis(:,[1 end]) = sds*segSis(:,[2 end-1]);

    segFes(2:end-1,2:end-1)   = segFes(2:end-1,2:end-1) + diff(segFes(:,2:end-1),2,1)./8 + diff(segFes(2:end-1,:),2,2)./8;
    segFes([1 end],:) = 0;
    segFes(:,[1 end]) = sds*segFes(:,[2 end-1]);

    segSil(2:end-1,2:end-1)   = segSil(2:end-1,2:end-1) + diff(segSil(:,2:end-1),2,1)./8 + diff(segSil(2:end-1,:),2,2)./8;
    segSil([1 end],:) = 0;
    segSil(:,[1 end]) = sds*segSil(:,[2 end-1]);

    segFel(2:end-1,2:end-1)   = segFel(2:end-1,2:end-1) + diff(segFel(:,2:end-1),2,1)./8 + diff(segFel(2:end-1,:),2,2)./8;
    segFel([1 end],:) = 0;
    segFel(:,[1 end]) = sds*segFel(:,[2 end-1]);

end

% update phase velocities
WlSi        = SOL.W + segSil;
UlSi        = SOL.U;
WsSi        = SOL.W + segSis;
UsSi        = SOL.U;
WlFe        = SOL.W + segFel;
UlFe        = SOL.U;
WsFe        = SOL.W + segFes;
UsFe        = SOL.U;

% update velocity divergence
Div_V(2:end-1,2:end-1) = ddz(SOL.W(:,2:end-1),NUM.h) ...                   % get velocity divergence
                       + ddx(SOL.U(2:end-1,:),NUM.h);
Div_V([1 end],:) = Div_V([2 end-1],:);                                     % apply boundary conditions
Div_V(:,[1 end]) = Div_V(:,[2 end-1]);

% update volume source
Div_rhoV =  + advection(MAT.rho.*CHM.xSi.*CHM.fSis,0.*SOL.U,segSis  ,NUM.h,NUM.h,NUM.ADVN,'flx') ...
            + advection(MAT.rho.*CHM.xSi.*CHM.fSil,0.*SOL.U,segSil  ,NUM.h,NUM.h,NUM.ADVN,'flx') ...
            + advection(MAT.rho.*CHM.xFe.*CHM.fFes,0.*SOL.U,segFes  ,NUM.h,NUM.h,NUM.ADVN,'flx') ...
            + advection(MAT.rho.*CHM.xFe.*CHM.fFel,0.*SOL.U,segFel  ,NUM.h,NUM.h,NUM.ADVN,'flx') ...
            + advection(MAT.rho                   ,   SOL.U,   SOL.W,NUM.h,NUM.h,NUM.ADVN,'flx');

VolSrc = -((MAT.rho-rhoo)/NUM.dt + (Div_rhoV - MAT.rho.*Div_V))./MAT.rho;

% set variable boundary conditions
UBG    = mean(mean(VolSrc(2:end-1,2:end-1)))./2 .* (NUM.L/2 - NUM.XU);
WBG    = mean(mean(VolSrc(2:end-1,2:end-1)))./2 .* (NUM.D/2 - NUM.ZW);

%% Calculate stress and strain rates
% update strain rates
DEF.exx(:,2:end-1) = diff(SOL.U,1,2)./NUM.h - Div_V(:,2:end-1)./3;                     % x-normal strain rate
DEF.exx([1 end],:) = DEF.exx([2 end-1],:);                                         % apply boundary conditions
DEF.exx(:,[1 end]) = DEF.exx(:,[2 end-1]);
DEF.ezz(2:end-1,:) = diff(SOL.W,1,1)./NUM.h - Div_V(2:end-1,:)./3;                        % z-normal strain rate
DEF.ezz([1 end],:) = DEF.ezz([2 end-1],:);                                         % apply boundary conditions
DEF.ezz(:,[1 end]) = DEF.ezz(:,[2 end-1]);
DEF.exz            = 1/2.*(diff(SOL.U,1,1)./NUM.h+diff(SOL.W,1,2)./NUM.h);                     % shear strain rate

% update stresses
DEF.txx = MAT.Eta .* DEF.exx;                                                          % x-normal stress
DEF.tzz = MAT.Eta .* DEF.ezz;                                                          % z-normal stress
DEF.txz = MAT.EtaC.* DEF.exz;   

% entropy production/heat dissipation rate
[grdTx,grdTz] = gradient(SOL.T,NUM.h);
EntProd = MAT.ks(2:end-1,2:end-1).*(grdTz(2:end-1,2:end-1).^2 + grdTx(2:end-1,2:end-1).^2)...
        + DEF.exx(2:end-1,2:end-1).*DEF.txx(2:end-1,2:end-1) + DEF.exx(2:end-1,2:end-1).*DEF.txx(2:end-1,2:end-1)...
        + 2.*(DEF.exz(1:end-1,1:end-1)+DEF.exz(2:end,1:end-1)+DEF.exz(1:end-1,2:end)+DEF.exz(2:end,2:end))./4 ...
           .*(DEF.txz(1:end-1,1:end-1)+DEF.txz(2:end,1:end-1)+DEF.txz(1:end-1,2:end)+DEF.txz(2:end,2:end))./4 ...
        +  CHM.xFe(2:end-1,2:end-1).*CHM.fFes(2:end-1,2:end-1).^2./Ksgr_x(2:end-1,2:end-1) .* ((segFes(1:end-1,2:end-1)+segFes(2:end,2:end-1))./2).^2  ...
        +  CHM.xFe(2:end-1,2:end-1).*CHM.fFel(2:end-1,2:end-1).^2./Ksgr_f(2:end-1,2:end-1) .* ((segFel(1:end-1,2:end-1)+segFel(2:end,2:end-1))./2).^2  ...
        +  CHM.xSi(2:end-1,2:end-1).*CHM.fSis(2:end-1,2:end-1).^2./Ksgr_x(2:end-1,2:end-1) .* ((segFes(1:end-1,2:end-1)+segFes(2:end,2:end-1))./2).^2 ;

%% update physical time step
dtadvn =  NUM.h/2   /max(abs([UlSi(:);WlSi(:);UsSi(:);WsSi(:);UlFe(:);WlFe(:);UsFe(:);WsFe(:)])); % stable timestep for advection
dtdiff = (NUM.h/2)^2/max(max(PHY.kTFe,PHY.kTSi)./MAT.rho(:)./PHY.Cp);                         % stable time step for T diffusion

NUM.dt = min(min(0.75*dtdiff,NUM.CFL * dtadvn),dtmax);                      % fraction of minimum stable time step

if NUM.dt==dtmax
    dtlimit = 'max step limited';
elseif dtdiff<dtadvn
    dtlimit = 'diffusion limited';
else
    dtlimit = 'advection limited';
end
  
