% print update header
% fprintf(1,'  ---  update materials & deformation \n');
tic;  % start clock on update

%% convert weight to volume fraction
% update T and chemical-dependent density
MAT.rhoSis = PHY.rhoSis.*(1 - PHY.aTSi.*(SOL.T-CHM.perTSi) - PHY.gCSi.*(CHM.csSi-CHM.cphsSi1));
MAT.rhoSil = PHY.rhoSil.*(1 - PHY.aTSi.*(SOL.T-CHM.perTSi) - PHY.gCSi.*(CHM.clSi-CHM.cphsSi1));
MAT.rhoFes = PHY.rhoFes.*(1 - PHY.aTFe.*(SOL.T-CHM.TFe1)   - PHY.gCFe.*(CHM.csFe-CHM.cphsFe1));
MAT.rhoFel = PHY.rhoFel.*(1 - PHY.aTFe.*(SOL.T-CHM.TFe1)   - PHY.gCFe.*(CHM.clFe-CHM.cphsFe1));

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

% get momentum permission
kv = permute(cat(3,etas,etaSil,etaFel),[3,1,2]);
Mv = permute(repmat(kv,1,1,1,3),[4,1,2,3])./permute(repmat(kv,1,1,1,3),[1,4,2,3]);
 
ff = max(1e-6,min(1-1e-6,permute(cat(3,MAT.phiSis+MAT.phiFes,MAT.phiSil,MAT.phiFel),[3,1,2])));
FF = permute(repmat(ff,1,1,1,3),[4,1,2,3]);
Sf = (FF./BBP).^(1./CCP);  Sf = Sf./sum(Sf,2);
Xf = sum(AAP.*Sf,2).*FF + (1-sum(AAP.*Sf,2)).*Sf;

thtv = squeeze(prod(Mv.^Xf,2));

% get momentum flux and transfer coefficients
Kv =    ff .*kv.*thtv;
Cv = (1-ff)./[PHY.dx;PHY.dm;PHY.df].^2.*Kv;

% compose effective viscosity, segregation coefficients
MAT.Eta  = squeeze(sum(Kv,1));                                             % effective magma viscosity
MAT.Eta  = (1/etamax + 1./(MAT.Eta + etamin)).^(-1);                       % limit viscosity range
MAT.EtaC = (MAT.Eta(1:end-1,1:end-1)+MAT.Eta(2:end,1:end-1) ...            % viscosity in cell corners
         +  MAT.Eta(1:end-1,2:end  )+MAT.Eta(2:end,2:end  ))./4;

Ksgr_x   = max(1e-15,min(1e-6,(MAT.phiSis+MAT.phiFes)./squeeze(Cv(1,:,:))));
Ksgr_m   = max(1e-15,min(1e-6, MAT.phiSil            ./squeeze(Cv(2,:,:))));
Ksgr_f   = max(1e-15,min(1e-6, MAT.phiFel            ./squeeze(Cv(3,:,:))));


%% other bulk thermochemical properties
% Update pressure
% SOL.Pt      = rhoRef.*PHY.gzP.*NUM.ZP + SOL.P;

% update heat capacities
% MAT.rhoCp = (MAT.rho.*CHM.xFe.*(CHM.fFes.*PHY.CpFes + CHM.fFel.*PHY.CpFel)...  % mixture sensible heat capacity density
%           +  MAT.rho.*CHM.xSi.*(CHM.fSis.*PHY.CpSis + CHM.fSil.*PHY.CpSil));
% MAT.rhoDs = (MAT.rho.*CHM.xFe.*CHM.fFel.*CHM.dEntrFe...                        % mixture latent heat capacity density
%           +  MAT.rho.*CHM.xSi.*CHM.fSil.*CHM.dEntrSi);
MAT.Ds    =  CHM.xFe.*CHM.fFel.*CHM.dEntrFe + CHM.xSi.*CHM.fSil.*CHM.dEntrSi;   % mixture entropy
MAT.kT    =  CHM.xFe.*PHY.kTFe + CHM.xSi.*PHY.kTSi;                            % magma thermal conductivity

%% calculate stokes settling velocity
sds = SOL.BCsides;      % side boundary type 
% also runge kutta advection time stepping
kappaseg = 100;
segSis      = ((MAT.rhoSis(1:end-1,:) + MAT.rhoSis(2:end,:))/2 ...
            -  (MAT.rho   (1:end-1,:) + MAT.rho   (2:end,:))/2)...
    .* (Ksgr_x    (1:end-1,:) .*Ksgr_x    (2:end,:)).^0.5  .*  PHY.gz; % silicate particle settling speed


segFes      = ((MAT.rhoFes(1:end-1,:) + MAT.rhoFes(2:end,:))/2 ...
            -  (MAT.rho   (1:end-1,:) + MAT.rho   (2:end,:))/2)...
    .* (Ksgr_x    (1:end-1,:) .*Ksgr_x    (2:end,:)).^0.5  .*  PHY.gz; % iron particle settling speed


segSil      = ((MAT.rhoSil(1:end-1,:) + MAT.rhoSil(2:end,:))/2 ...
            -  (MAT.rho   (1:end-1,:) + MAT.rho   (2:end,:))/2)...
    .* (Ksgr_m    (1:end-1,:) .*Ksgr_m    (2:end,:)).^0.5  .*  PHY.gz.*Siseg; % silicate particle settling speed


segFel      = ((MAT.rhoFel(1:end-1,:) + MAT.rhoFel(2:end,:))/2 ...
            -  (MAT.rho   (1:end-1,:) + MAT.rho   (2:end,:))/2)...
    .* (Ksgr_f    (1:end-1,:) .*Ksgr_f    (2:end,:)).^0.5  .*  PHY.gz; % iron particle settling speed


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


%% update physical time step
dtadvn =  NUM.h/2   /max(abs([UlSi(:);WlSi(:);UsSi(:);WsSi(:);UlFe(:);WlFe(:);UsFe(:);WsFe(:)])); % stable timestep for advection
dtdiff = (NUM.h/2)^2/max(MAT.kT(:)./MAT.rho(:)./PHY.Cp);                         % stable time step for T diffusion

NUM.dt = min(min(0.75*dtdiff,NUM.CFL * dtadvn),dtmax);                      % fraction of minimum stable time step

if NUM.dt==dtmax
    dtlimit = 'max step limited';
elseif dtdiff<dtadvn
    dtlimit = 'diffusion limited';
else
    dtlimit = 'advection limited';
end
  
