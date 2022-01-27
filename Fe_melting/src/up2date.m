% print update header
fprintf(1,'  ---  update materials & deformation \n');
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
% MAT.EtalSi = zeros(size(SOL.phiSil)) + PHY.EtalSi0;
% MAT.EtasSi = zeros(size(SOL.phiSis)) + PHY.EtasSi0;
% MAT.EtasFe = zeros(size(SOL.phiFes)) + PHY.EtasFe0;
% % as a function of permission weights
% kv = permute(cat(3,PHY.Eta0,PHY.EtasSi0,PHY.EtasFe0),[3,1,2]);
% Mv = permute(repmat(kv,1,1,1,3),[4,1,2,3])./permute(repmat(kv,1,1,1,3),[1,4,2,3]);

% legacy
MAT.EtaS                = (1-SOL.phiFes./0.5).^(-2)...  % iron crystal suspension
                        .*(1-SOL.phiSis./0.5).^(-2);    % silicate crystal suspension
% MAT.Mu                  = -4.13 +3703./(SOL.T-761.7);   % T-dependent melt viscosity
MAT.Eta                 = MAT.Eta0.*MAT.EtaS;           % mixture viscosity
MT.ETA(SOL.phiSil<= 0.5)= 1e17;
% MAT.Eta                 = MAT.Mu.*MAT.EtaS;           % mixture viscosity

MAT.EtaC                = (MAT.Eta(1:end-1,1:end-1) ...
                        +  MAT.Eta(2:end  ,1:end-1) ...
                        +  MAT.Eta(1:end-1,2:end  ) ...
                        +  MAT.Eta(2:end  ,2:end  ))/4;                     % interpolate to corner nodes

% %% update T and chemical-dependent density
% MAT.rhoUt               = (MAT.rhot(:,1:end-1)+MAT.rhot(:,2:end))./2;     % density on the vx nodes
% MAT.rhoWt               = (MAT.rhot(1:end-1,:)+MAT.rhot(2:end,:))./2;     % density on the vz nodes                 
%                     
% % silicate
% MAT.RhoSi               = PHY.RhoSi0.* (1 - MAT.aT.*(SOL.T-SOL.T0));        % density on centre nodes
% MAT.RhoUSi              = (MAT.RhoSi(:,1:end-1)+MAT.RhoSi(:,2:end))./2;   % density on the vx nodes
% MAT.RhoWSi              = (MAT.RhoSi(1:end-1,:)+MAT.RhoSi(2:end,:))./2;   % density on the vz nodes
% 
% % iron droplets
% MAT.RhoFe               = PHY.RhoFe0.*(1 - MAT.aT.*(SOL.T-SOL.T0));        % density on centre nodes
% MAT.RhoUFe              = (MAT.RhoFe(:,1:end-1)+MAT.RhoFe(:,2:end))./2;   % density on the vx nodes
% MAT.RhoWFe              = (MAT.RhoFe(1:end-1,:)+MAT.RhoFe(2:end,:))./2;   % density on the vz nodes
% 
% % density contrast (between Fe and bulk)
% MAT.dRho                = MAT.RhoFe - MAT.Rhot;

%% other bulk thermochemical properties
% Update pressure
% rhoRef  = mean(mean(MAT.rhot(2:end-1,2:end-1)));
SOL.Pt  = rhoRef.*PHY.gzP.*NUM.ZP + SOL.P;

MAT.rhoCpt  = (MAT.rhot.*CHM.xFe.*(CHM.fFes.*PHY.CpFes + CHM.fFel.*PHY.CpFel)...
            +  MAT.rhot.*CHM.xSi.*(CHM.fSis.*PHY.CpSis + CHM.fSil.*PHY.CpSil));
MAT.kT      = PHY.kTFe.*(SOL.phiFes + SOL.phiFel)...
            + PHY.kTSi.*(SOL.phiSis + SOL.phiSil);

%% calculate stokes settling velocity  
if SOL.BCSeg==2; sds = -1;      % no slip
else;            sds = +1; end  % free slip


segSis      = 2/9 .* ((MAT.rhoSis(1:end-1,:)+ MAT.rhoSis(2:end,:))./2 ...
                     -(MAT.rhot(1:end-1,:)  + MAT.rhot(2:end,:))/2)...
                  .*  PHY.gz*PHY.d^2 ...
                  ./ ((MAT.Eta(1:end-1,:)   + MAT.Eta(2:end,:) )/2); % crystal settling speed
segSis([1 end],:) = 0;
segSis(:,[1 end]) = sds*segSis(:,[2 end-1]);

segFes      = 2/9 .* ((MAT.rhoFes(1:end-1,:)+ MAT.rhoFes(2:end,:))./2      ...
                     -(MAT.rhot(1:end-1,:)  + MAT.rhot(2:end,:))/2)...
                  .*  PHY.gz*PHY.d^2  ...
                  ./ ((MAT.Eta(1:end-1,:)   + MAT.Eta(2:end,:) )/2); % crystal settling speed
segFes([1 end],:) = 0;
segFes(:,[1 end]) = sds*segFes(:,[2 end-1]);

segFel      = 2/9 .* ((MAT.rhoFel(1:end-1,:)+ MAT.rhoFel(2:end,:))./2 ...
                     -(MAT.rhot(1:end-1,:)  + MAT.rhot(2:end,:))/2)...
                  .*  PHY.gz*PHY.d^2 ...
                  ./((MAT.Eta(1:end-1,:)   + MAT.Eta(2:end,:) )/2); % crystal settling speed
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
VolSrc = -((MAT.rhot-MAT.rhoo)./NUM.dt + (Div_rhov - MAT.rhot.*Div_V))./MAT.rhot;

dVoldt = mean(mean(VolSrc(2:end-1,2:end-1)));
VolSrc = VolSrc - dVoldt;


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

NUM.dt              = NUM.CFL * min(dtdiff,dtadvn);                     % fraction of minimum stable time step
if dtdiff<dtadvn
    disp('diffusion regime')
else
    disp('advection regime')
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