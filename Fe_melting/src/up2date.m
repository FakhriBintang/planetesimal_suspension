% print update header
% fprintf(1,'  ---  update materials & deformation \n');
tic;  % start clock on update

%% convert weight to volume fraction
%% update T and chemical-dependent density
MAT.rhoSis              = PHY.rhoSis.*(1 - MAT.aTSi.*(SOL.T-SOL.T0) - PHY.gammaSi.*(CHM.csSi-CHM.cphsSi1));  
MAT.rhoSil              = PHY.rhoSil.*(1 - MAT.aTSi.*(SOL.T-SOL.T0) - PHY.gammaSi.*(CHM.clSi-CHM.cphsSi1));
MAT.rhoFes              = PHY.rhoFes.*(1 - MAT.aTFe.*(SOL.T-SOL.T0) - PHY.gammaFe.*(CHM.csFe-CHM.cphsFe1));  
MAT.rhoFel              = PHY.rhoFel.*(1 - MAT.aTFe.*(SOL.T-SOL.T0) - PHY.gammaFe.*(CHM.clFe-CHM.cphsFe1));

MAT.rhot   = 1./ (CHM.xFe.*CHM.fFes./MAT.rhoFes + CHM.xFe.*CHM.fFel./MAT.rhoFel + CHM.xSi.*CHM.fSis./MAT.rhoSis + CHM.xSi.*CHM.fSil./MAT.rhoSil);

SOL.phiFes     = CHM.xFe.* CHM.fFes .* MAT.rhot ./ MAT.rhoFes;
SOL.phiFel     = CHM.xFe.* CHM.fFel .* MAT.rhot ./ MAT.rhoFel;
SOL.phiSis     = CHM.xSi.* CHM.fSis .* MAT.rhot ./ MAT.rhoSis;
SOL.phiSil     = CHM.xSi.* CHM.fSil .* MAT.rhot ./ MAT.rhoSil; 

% MAT.rhoSis              = PHY.rhoSis.*(1 - MAT.aTSi.*(SOL.T-SOL.T0) - PHY.gammaSi.*cSi);  
% MAT.rhoSil              = PHY.rhoSil.*(1 - MAT.aTSi.*(SOL.T-SOL.T0) - PHY.gammaSi.*cSi);
% MAT.rhoFes              = PHY.rhoFes.*(1 - MAT.aTFe.*(SOL.T-SOL.T0) - PHY.gammaFe.*cFe);  
% MAT.rhoFel              = PHY.rhoFel.*(1 - MAT.aTFe.*(SOL.T-SOL.T0) - PHY.gammaFe.*cFe);
% for loop = 1:5
% % CHM.fFel = 1-CHM.fFes; CHM.fSil = 1-CHM.fSis;
% 
% SOL.phiFes = max(TINY,min(1-TINY, CHM.xFe.*CHM.fFes.*MAT.rhot./MAT.rhoFes ))./(SOL.phiFes+SOL.phiFel+SOL.phiSis+SOL.phiSil);
% SOL.phiFel = max(TINY,min(1-TINY, CHM.xFe.*CHM.fFel.*MAT.rhot./MAT.rhoFel ))./(SOL.phiFes+SOL.phiFel+SOL.phiSis+SOL.phiSil);
% SOL.phiSis = max(TINY,min(1-TINY, CHM.xSi.*CHM.fSis.*MAT.rhot./MAT.rhoSis ))./(SOL.phiFes+SOL.phiFel+SOL.phiSis+SOL.phiSil);
% SOL.phiSil = max(TINY,min(1-TINY, CHM.xSi.*CHM.fSil.*MAT.rhot./MAT.rhoSil ))./(SOL.phiFes+SOL.phiFel+SOL.phiSis+SOL.phiSil);
% 
% MAT.rhot                = SOL.phiFes.*MAT.rhoFes...
%                         + SOL.phiFel.*MAT.rhoFel...
%                         + SOL.phiSis.*MAT.rhoSis...
%                         + SOL.phiSil.*MAT.rhoSil;
% % figure(200); 
% % subplot(2,1,1); imagesc(SOL.phiFel+SOL.phiFes+SOL.phiSil+SOL.phiSis); colorbar; title('phisum')
% % subplot(2,1,2); imagesc(MAT.rhot); colorbar
% % drawnow
% % CHM.XFe = MAT.rhot.*CHM.xFe; CHM.XSi = MAT.rhot.*CHM.xSi;        
% end

%% update viscosity
etas = zeros(size(SOL.phiSis)) + PHY.Etasol0; 
etam = PHY.Eta0 .* exp(Em./(8.3145.*(SOL.T+273.15))-Em./(8.3145.*((CHM.TFe1+CHM.TFe2)/2+273.15)));
etaFe= zeros(size(SOL.phiFel)) + PHY.EtalFe0;
% get permission weights
kv = permute(cat(3,etas,etam,etaFe),[3,1,2]);
Mv = permute(repmat(kv,1,1,1,3),[4,1,2,3])./permute(repmat(kv,1,1,1,3),[1,4,2,3]);
 
ff = max(1e-16,permute(cat(3,SOL.phiSis+SOL.phiFes,SOL.phiSil,SOL.phiFel),[3,1,2]));
FF = permute(repmat(ff,1,1,1,3),[4,1,2,3]);
Sf = (FF./BBP).^(1./CCP);  Sf = Sf./sum(Sf,2);
Xf = sum(AAP.*Sf,2).*FF + (1-sum(AAP.*Sf,2)).*Sf;

% get momentum and volume permissions
thtv = squeeze(prod(Mv.^Xf,2));

% get momentum and volume flux and transfer coefficients
Kv =    ff .*kv.*thtv;
Cv = (1-ff)./[PHY.dx;PHY.dm;PHY.df].^2.*Kv;

% compose effective viscosity, segregation coefficients
Kvb = squeeze(sum(Kv,1));                                                  % effective magma viscosity
if NUM.step>0; MAT.Eta = (Kvb + Kvbo)/2; 
else;      MAT.Eta =  Kvb;  end
MAT.Eta   = max(etamin,min(etamax,MAT.Eta));                                       % limit viscosity range
MAT.EtaC  = (MAT.Eta(1:end-1,1:end-1)+MAT.Eta(2:end,1:end-1) ...                       % viscosity in cell corners
          +  MAT.Eta(1:end-1,2:end  )+MAT.Eta(2:end,2:end  ))./4;

Ksgr_x = max(1e-18,min(1e-6,SOL.phiSis+SOL.phiFes./squeeze(Cv(1,:,:))));
Ksgr_m = max(1e-18,min(1e-6,SOL.phiSil ./squeeze(Cv(2,:,:))));
Ksgr_f = max(1e-18,min(1e-6,SOL.phiFel ./squeeze(Cv(3,:,:))));

%% other bulk thermochemical properties
% Update pressure
% rhoRef  = mean(mean(MAT.rhot(2:end-1,2:end-1)));
SOL.Pt      = rhoRef.*PHY.gzP.*NUM.ZP + SOL.P;

MAT.rhoCpt  = (MAT.rhot.*CHM.xFe.*(CHM.fFes.*PHY.CpFes + CHM.fFel.*PHY.CpFel)...
            +  MAT.rhot.*CHM.xSi.*(CHM.fSis.*PHY.CpSis + CHM.fSil.*PHY.CpSil));
MAT.kT      = PHY.kTFe.*(SOL.phiFes + SOL.phiFel)...
            + PHY.kTSi.*(SOL.phiSis + SOL.phiSil);

%% calculate stokes settling velocity  
if SOL.BCSeg==2; sds = -1;      % no slip
else;            sds = +1; end  % free slip


segSis      = ((MAT.rhoSis(1:end-1,:)+ MAT.rhoSis(2:end,:))./2 ...
              -(MAT.rhot(1:end-1,:)  + MAT.rhot(2:end,:))  ./2)...
            .*((Ksgr_x(1:end-1,:)    +Ksgr_x(2:end,:))/2)  .*  PHY.gz; % crystal settling speed
segSis([1 end],:) = 0;
segSis(:,[1 end]) = sds*segSis(:,[2 end-1]);

segFes      = ((MAT.rhoFes(1:end-1,:)+ MAT.rhoFes(2:end,:))./2 ...
              -(MAT.rhot(1:end-1,:)  + MAT.rhot(2:end,:))  ./2)...
            .*((Ksgr_x(1:end-1,:)    +Ksgr_x(2:end,:))/2)  .*  PHY.gz; % Fe crystal settling speed
segFes([1 end],:) = 0;
segFes(:,[1 end]) = sds*segFes(:,[2 end-1]);

segFel      = ((MAT.rhoFel(1:end-1,:)+ MAT.rhoFel(2:end,:))./2 ...
              -(MAT.rhot(1:end-1,:)  + MAT.rhot(2:end,:))  ./2)...
            .*((Ksgr_f(1:end-1,:)    +Ksgr_f(2:end,:))/2)  .*  PHY.gz; % droplet settling speed
segFel([1 end],:) = 0;
segFel(:,[1 end]) = sds*segFel(:,[2 end-1]);

segSil      = ((MAT.rhoSil(1:end-1,:)+ MAT.rhoSil(2:end,:))./2 ...
              -(MAT.rhot(1:end-1,:)  + MAT.rhot(2:end,:))  ./2)...
            .*((Ksgr_m(1:end-1,:)    +Ksgr_m(2:end,:))/2)  .*  PHY.gz; % droplet settling speed
segSil([1 end],:) = 0;
segSil(:,[1 end]) = sds*segSil(:,[2 end-1]);

% update phase velocities
WlSi        = SOL.W;
UlSi        = SOL.U;
WsSi        = SOL.W + segSis;
UsSi        = SOL.U;
WlFe        = SOL.W + segFel;
UlFe        = SOL.U;
WsFe        = SOL.W + segFes;
UsFe        = SOL.U;
% if ~mod(NUM.step,round(2*RUN.nup/NUM.CFL))                          % only update deformation when fluid mechanics solved

% end

% update velocity divergence
Div_V(2:end-1,2:end-1) = ddz(SOL.W(:,2:end-1),NUM.h) ...                   % get velocity divergence
                       + ddx(SOL.U(2:end-1,:),NUM.h);
Div_V([1 end],:) = Div_V([2 end-1],:);                                     % apply boundary conditions
Div_V(:,[1 end]) = Div_V(:,[2 end-1]);

% update volume source
Div_rhov =  + advection(MAT.rhot.*CHM.xSi.*CHM.fSis,0.*SOL.U,segSis  ,NUM.h,NUM.h,ADVN,'flx') ...
            + advection(MAT.rhot.*CHM.xSi.*CHM.fSil,0.*SOL.U,0.*SOL.W,NUM.h,NUM.h,ADVN,'flx') ...
            + advection(MAT.rhot.*CHM.xFe.*CHM.fFes,0.*SOL.U,segFes  ,NUM.h,NUM.h,ADVN,'flx') ...
            + advection(MAT.rhot.*CHM.xFe.*CHM.fFel,0.*SOL.U,segFel  ,NUM.h,NUM.h,ADVN,'flx') ...
            + advection(MAT.rhot           ,   SOL.U,   SOL.W,NUM.h,NUM.h,ADVN,'flx');

% VolSrc = -((MAT.rhot-MAT.rhoo)./NUM.dt + (Div_rhov - MAT.rhot.*Div_V + Div_rhoVo)/2)./(MAT.rhot/2);
VolSrc  = -((MAT.rhot-MAT.rhoo)./NUM.dt + (Div_rhov - MAT.rhot.*Div_V))./MAT.rhot;
% set variable boundary conditions
UBG     = mean(mean(VolSrc(2:end-1,2:end-1)))./2 .* (NUM.L/2 - NUM.XU);
WBG     = mean(mean(VolSrc(2:end-1,2:end-1)))./2 .* (NUM.D/2 - NUM.ZW);

% legacy
% dVoldt = mean(mean(VolSrc(2:end-1,2:end-1)));
% VolSrc = VolSrc - dVoldt;


%% update physical time step
VR                  = abs([SOL.U(:);SOL.W(:)]);                         % reference velocity
dtadvR              = NUM.h/2/max(VR);
VSis                = abs([SOL.U(:);SOL.W(:)+segSis(:)]);              % silicate crystal velocity
dtadvSis            = NUM.h/2/max(VSis);                                  %
VFes                = abs([SOL.U(:);SOL.W(:)+segFes(:)]);              % Fe solid velocity
dtadvFes            = NUM.h/2/max(VFes);                                  %
VFel                = abs([SOL.U(:);SOL.W(:)+segFel(:)]);              % Fe solid velocity
dtadvFel            = NUM.h/2/max(VFel);                                  %

dtadvn              = min([dtadvR dtadvSis dtadvFes dtadvFel]);                               % Stable timestep for advection

kappa               = MAT.kT./rhoCpt;
kappa               = kappa(:);
dtdiff              = (NUM.h/2)^2 / max(kappa);            % stable time step for T diffusion

NUM.dt              = min(NUM.CFL * min(dtdiff,dtadvn),dtmax);                     % fraction of minimum stable time step
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