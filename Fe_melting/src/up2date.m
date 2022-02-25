% print update header
% fprintf(1,'  ---  update materials & deformation \n');
tic;  % start clock on update

%% convert weight to volume fraction
% update T and chemical-dependent density
MAT.rhoSis = PHY.rhoSis.*(1 - PHY.aTSi.*(SOL.T-SOL.T0) - PHY.gCSi.*(CHM.csSi-CHM.cphsSi1));
MAT.rhoSil = PHY.rhoSil.*(1 - PHY.aTSi.*(SOL.T-SOL.T0) - PHY.gCSi.*(CHM.clSi-CHM.cphsSi1));
MAT.rhoFes = PHY.rhoFes.*(1 - PHY.aTFe.*(SOL.T-SOL.T0) - PHY.gCFe.*(CHM.csFe-CHM.cphsFe1));
MAT.rhoFel = PHY.rhoFel.*(1 - PHY.aTFe.*(SOL.T-SOL.T0) - PHY.gCFe.*(CHM.clFe-CHM.cphsFe1));

% update mixture density
MAT.rho    = 1./(CHM.xFe.*CHM.fFes./MAT.rhoFes + CHM.xFe.*CHM.fFel./MAT.rhoFel + CHM.xSi.*CHM.fSis./MAT.rhoSis + CHM.xSi.*CHM.fSil./MAT.rhoSil);

% update volume fractions
SOL.phiFes = CHM.xFe.* CHM.fFes .* MAT.rho ./ MAT.rhoFes;
SOL.phiFel = CHM.xFe.* CHM.fFel .* MAT.rho ./ MAT.rhoFel;
SOL.phiSis = CHM.xSi.* CHM.fSis .* MAT.rho ./ MAT.rhoSis;
SOL.phiSil = CHM.xSi.* CHM.fSil .* MAT.rho ./ MAT.rhoSil; 


%% update viscosity
% get pure phase densities
etas   = zeros(size(SOL.phiSis)) + PHY.EtaSol0; 
etaSil = PHY.EtaSil0 .* exp(Em./(8.3145.*(SOL.T+273.15))-Em./(8.3145.*((CHM.TFe1+CHM.TFe2)/2+273.15)));
etaFel = zeros(size(SOL.phiFel)) + PHY.EtaFel0;

% get momentum permission
kv = permute(cat(3,etas,etaSil,etaFel),[3,1,2]);
Mv = permute(repmat(kv,1,1,1,3),[4,1,2,3])./permute(repmat(kv,1,1,1,3),[1,4,2,3]);
 
ff = max(1e-16,permute(cat(3,SOL.phiSis+SOL.phiFes,SOL.phiSil,SOL.phiFel),[3,1,2]));
FF = permute(repmat(ff,1,1,1,3),[4,1,2,3]);
Sf = (FF./BBP).^(1./CCP);  Sf = Sf./sum(Sf,2);
Xf = sum(AAP.*Sf,2).*FF + (1-sum(AAP.*Sf,2)).*Sf;

thtv = squeeze(prod(Mv.^Xf,2));

% get momentum flux and transfer coefficients
Kv =    ff .*kv.*thtv;
Cv = (1-ff)./[PHY.dx;PHY.dm;PHY.df].^2.*Kv;

% compose effective viscosity, segregation coefficients
MAT.Eta  = squeeze(sum(Kv,1));                                             % effective magma viscosity
MAT.Eta  = max(etamin,min(etamax,MAT.Eta));                                % limit viscosity range
MAT.EtaC = (MAT.Eta(1:end-1,1:end-1)+MAT.Eta(2:end,1:end-1) ...            % viscosity in cell corners
         +  MAT.Eta(1:end-1,2:end  )+MAT.Eta(2:end,2:end  ))./4;

Ksgr_x   = max(1e-18,min(1e-6,(SOL.phiSis+SOL.phiFes)./squeeze(Cv(1,:,:))));
Ksgr_m   = max(1e-18,min(1e-6, SOL.phiSil            ./squeeze(Cv(2,:,:))));
Ksgr_f   = max(1e-18,min(1e-6, SOL.phiFel            ./squeeze(Cv(3,:,:))));


%% other bulk thermochemical properties
% Update pressure
% SOL.Pt      = rhoRef.*PHY.gzP.*NUM.ZP + SOL.P;

% update heat capacities
MAT.rhoCp = (MAT.rho.*CHM.xFe.*(CHM.fFes.*PHY.CpFes + CHM.fFel.*PHY.CpFel)...  % mixture sensible heat capacity density
          +  MAT.rho.*CHM.xSi.*(CHM.fSis.*PHY.CpSis + CHM.fSil.*PHY.CpSil));
MAT.rhoDs = (MAT.rho.*CHM.xFe.*CHM.fFel.*CHM.dEntrFe...                        % mixture latent heat capacity density
          +  MAT.rho.*CHM.xSi.*CHM.fSil.*CHM.dEntrSi);
MAT.kT    =  CHM.xFe.*PHY.kTFe + CHM.xSi.*PHY.kTSi;                            % magma thermal conductivity

        
%% calculate stokes settling velocity  
sds = SOL.BCsides;      % side boundary type

segSis      = ((MAT.rhoSis(1:end-1,:) + MAT.rhoSis(2:end,:))/2 ...
            -  (MAT.rho   (1:end-1,:) + MAT.rho   (2:end,:))/2)...
            .*((Ksgr_x    (1:end-1,:) + Ksgr_x    (2:end,:))/2)  .*  PHY.gz; % silicate particle settling speed
segSis([1 end],:) = 0;
segSis(:,[1 end]) = sds*segSis(:,[2 end-1]);

segFes      = ((MAT.rhoFes(1:end-1,:) + MAT.rhoFes(2:end,:))/2 ...
            -  (MAT.rho   (1:end-1,:) + MAT.rho   (2:end,:))/2)...
            .*((Ksgr_x    (1:end-1,:) + Ksgr_x    (2:end,:))/2)  .*  PHY.gz; % iron particle settling speed
segFes([1 end],:) = 0;
segFes(:,[1 end]) = sds*segFes(:,[2 end-1]);

segSil      = ((MAT.rhoSil(1:end-1,:) + MAT.rhoSil(2:end,:))/2 ...
            -  (MAT.rho   (1:end-1,:) + MAT.rho   (2:end,:))/2)...
            .*((Ksgr_m    (1:end-1,:) + Ksgr_m    (2:end,:))/2)  .*  PHY.gz; % silicate particle settling speed
segSil([1 end],:) = 0;
segSil(:,[1 end]) = sds*segSil(:,[2 end-1]);

segFel      = ((MAT.rhoFel(1:end-1,:) + MAT.rhoFel(2:end,:))/2 ...
            -  (MAT.rho   (1:end-1,:) + MAT.rho   (2:end,:))/2)...
            .*((Ksgr_f    (1:end-1,:) + Ksgr_f    (2:end,:))/2)  .*  PHY.gz; % iron particle settling speed
segFel([1 end],:) = 0;
segFel(:,[1 end]) = sds*segFel(:,[2 end-1]);

% update phase velocities
WlSi        = SOL.W;
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

VolSrc = -((MAT.rho-rhoo)/NUM.dt + (Div_rhoV - MAT.rho.*Div_V + Div_rhoVo)/2)./MAT.rho;
% VolSrc = -((MAT.rho-rhoo)./NUM.dt + (Div_rhoV - MAT.rho.*Div_V))./MAT.rho;

% set variable boundary conditions
UBG    = mean(mean(VolSrc(2:end-1,2:end-1)))./2 .* (NUM.L/2 - NUM.XU);
WBG    = mean(mean(VolSrc(2:end-1,2:end-1)))./2 .* (NUM.D/2 - NUM.ZW);


%% update physical time step
dtadvn =  NUM.h/2   /max(abs([UlSi(:);WlSi(:);UsSi(:);WsSi(:);UlFe(:);WlFe(:);UsFe(:);WsFe(:)])); % stable timestep for advection
dtdiff = (NUM.h/2)^2/max(MAT.kT(:)./MAT.rhoCp(:));                         % stable time step for T diffusion

NUM.dt = min(NUM.CFL * min(dtdiff,dtadvn),dtmax);                          % fraction of minimum stable time step
if dtdiff<dtadvn
    dtlimit = 'diffusion limited';
else
    dtlimit = 'advection limited';
end
    
%     %% update strain-rate components
%     % get volumetric strain-rate (velocity divergence)
%     DEF.ups(2:end-1,2:end-1) = diff(SOL.U(2:end-1,:),1,2)./NUM.h ...
%                              + diff(SOL.W(:,2:end-1),1,1)./NUM.h;           % velocity divergence
%     DEF.ups([1 end],:)       = DEF.ups([2 end-1],:);                        % apply boundary conditions
%     DEF.ups(:,[1 end])       = DEF.ups(:,[2 end-1]);
%     
%     % get deviatoric strain rates
%     DEF.exx(:,2:end-1)      = diff(SOL.U,1,2)./NUM.h ...
%                             - DEF.ups(:,2:end-1)./3;                        % x-normal strain rate
%     DEF.exx([1 end],:)      = DEF.exx([2 end-1],:);                         % apply boundary conditions
%     DEF.exx(:,[1 end])      = DEF.exx(:,[2 end-1]);
%     
%     DEF.ezz(2:end-1,:)      = diff(SOL.W,1,1)./NUM.h ...
%                             - DEF.ups(2:end-1,:)./3;                        % z-normal strain rate
%     DEF.ezz([1 end],:)      =  DEF.ezz([2 end-1],:);                        % apply boundary conditions
%     DEF.ezz(:,[1 end])      =  DEF.ezz(:,[2 end-1]);
%     
%     DEF.exz                 = (diff(SOL.U,1,1)./NUM.h ...
%                             +  diff(SOL.W,1,2)./NUM.h)/2;                   % shear strain rate
%     
%     % update strain-rate magnitude
%     DEF.eII(2:end-1,2:end-1)= (  (DEF.exx(2:end-1,2:end-1).^2 + DEF.ezz(2:end-1,2:end-1).^2 ...
%                             + 2.*(DEF.exz(1:end-1,1:end-1).^2 + DEF.exz(2:end,1:end-1).^2 ...
%                             +     DEF.exz(1:end-1,2:end  ).^2 + DEF.exz(2:end,2:end  ).^2)/4 )/2).^0.5 + 1e-16;
%     DEF.eII([1 end],:)      =     DEF.eII([2 end-1],:);                   	% apply boundaries
%     DEF.eII(:,[1 end])      =     DEF.eII(:,[2 end-1]);
%     
%     
%     %% update stress components
%     DEF.txx = MAT.Eta  .* DEF.exx;                                          % x-normal stress
%     DEF.tzz = MAT.Eta  .* DEF.ezz;                                          % z-normal stress
%     DEF.txz = MAT.EtaC .* DEF.exz;                                          % xz-shear stress
%     
%     % update strain-rate magnitude
%     DEF.tII(2:end-1,2:end-1)= (  (DEF.txx(2:end-1,2:end-1).^2 + DEF.tzz(2:end-1,2:end-1).^2 ...
%                             + 2.*(DEF.txz(1:end-1,1:end-1).^2 + DEF.txz(2:end,1:end-1).^2 ...
%                             +     DEF.txz(1:end-1,2:end  ).^2 + DEF.txz(2:end,2:end  ).^2)/4 )/2).^0.5 + 1e-16;
%     DEF.tII([1 end],:)      =     DEF.tII([2 end-1],:);                     % apply boundaries
%     DEF.tII(:,[1 end])      =     DEF.tII(:,[2 end-1]);
% end
% %% update heat source fields
%     % update shear heating
%     SOL.Hs = 2.*DEF.eII.*DEF.tII;
%     
%     % update adiabatic heating
%     SOL.Ha = ((SOL.WP.*PHY.gzP + SOL.UP.*PHY.gxP) .* MAT.Rhot + SOL.seg.*PHY.gzP).*MAT.aT.*SOL.T;