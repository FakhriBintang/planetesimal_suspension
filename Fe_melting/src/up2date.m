% print update header
fprintf(1,'  ---  update materials & deformation \n');
tic;  % start clock on update

%% convert weight to volume fraction
%% update T and chemical-dependent density
MAT.rhoSis              = PHY.rhoSis.*(1 - MAT.aTSi.*(SOL.T-SOL.T0) - PHY.gammaSi.*cSi);  
MAT.rhoSil              = PHY.rhoSil.*(1 - MAT.aTSi.*(SOL.T-SOL.T0) - PHY.gammaSi.*cSi);
MAT.rhoFes              = PHY.rhoFes.*(1 - MAT.aTFe.*(SOL.T-SOL.T0) - PHY.gammaFe.*cFe);  
MAT.rhoFel              = PHY.rhoFel.*(1 - MAT.aTFe.*(SOL.T-SOL.T0) - PHY.gammaFe.*cFe);

MAT.rhot   = 1./ (xFe.*fFes./MAT.rhoFes + xFe.*fFel./MAT.rhoFel + xSi.*fSis./MAT.rhoSis + xSi.*fSil./MAT.rhoSil);

phiFes     = xFe.* fFes .* MAT.rhot ./ MAT.rhoFes;
phiFel     = xFe.* fFel .* MAT.rhot ./ MAT.rhoFel;
phiSis     = xSi.* fSis .* MAT.rhot ./ MAT.rhoSis;
phiSil     = xSi.* fSil .* MAT.rhot ./ MAT.rhoSil; 

% MAT.rhoSis              = PHY.rhoSis.*(1 - MAT.aTSi.*(SOL.T-SOL.T0) - PHY.gammaSi.*cSi);  
% MAT.rhoSil              = PHY.rhoSil.*(1 - MAT.aTSi.*(SOL.T-SOL.T0) - PHY.gammaSi.*cSi);
% MAT.rhoFes              = PHY.rhoFes.*(1 - MAT.aTFe.*(SOL.T-SOL.T0) - PHY.gammaFe.*cFe);  
% MAT.rhoFel              = PHY.rhoFel.*(1 - MAT.aTFe.*(SOL.T-SOL.T0) - PHY.gammaFe.*cFe);
% for loop = 1:5
% % fFel = 1-fFes; fSil = 1-fSis;
% 
% phiFes = max(TINY,min(1-TINY, xFe.*fFes.*MAT.rhot./MAT.rhoFes ))./(phiFes+phiFel+phiSis+phiSil);
% phiFel = max(TINY,min(1-TINY, xFe.*fFel.*MAT.rhot./MAT.rhoFel ))./(phiFes+phiFel+phiSis+phiSil);
% phiSis = max(TINY,min(1-TINY, xSi.*fSis.*MAT.rhot./MAT.rhoSis ))./(phiFes+phiFel+phiSis+phiSil);
% phiSil = max(TINY,min(1-TINY, xSi.*fSil.*MAT.rhot./MAT.rhoSil ))./(phiFes+phiFel+phiSis+phiSil);
% 
% MAT.rhot                = phiFes.*MAT.rhoFes...
%                         + phiFel.*MAT.rhoFel...
%                         + phiSis.*MAT.rhoSis...
%                         + phiSil.*MAT.rhoSil;
% % figure(200); 
% % subplot(2,1,1); imagesc(phiFel+phiFes+phiSil+phiSis); colorbar; title('phisum')
% % subplot(2,1,2); imagesc(MAT.rhot); colorbar
% % drawnow
% % XFe = MAT.rhot.*xFe; XSi = MAT.rhot.*xSi;        
% end

% update viscosity
MAT.EtaS                = (1-phiFes./0.5).^(-2)...  % iron crystal suspension
                        .*(1-phiSis./0.5).^(-2);    % silicate crystal suspension
% MAT.Mu                  = -4.13 +3703./(SOL.T-761.7);   % T-dependent melt viscosity
MAT.Eta                 = MAT.Eta0.*MAT.EtaS;           % mixture viscosity
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
SOL.Pt = MAT.rhot.*PHY.gzP.*NUM.ZP + SOL.P;

rhoRef  = mean(mean(MAT.rhot(2:end-1,2:end-1)));
MAT.rhoCpt  = (XFe.*(fFes.*PHY.CpFes + fFel.*PHY.CpFel)...
            +  XSi.*(fSis.*PHY.CpSis + fSil.*PHY.CpSil));
MAT.kT      = PHY.kTFe.*(phiFes + phiFel)...
            + PHY.kTSi.*(phiSis + phiSil);

%% calculate stokes settling velocity  
if SOL.BCSeg==2; sds = -1;      % no slip
else;            sds = +1; end  % free slip


segSis      = -2/9 .* ((MAT.rhoSis(1:end-1,:)+ MAT.rhoSis(2:end,:))./2 ...
                     -(MAT.rhot(1:end-1,:)  + MAT.rhot(2:end,:))/2)...
                  .*  PHY.gz*PHY.d^2.* (phiSis(1:end-1,:)+phiSis(2:end,:))./2 ...
                  ./ ((MAT.Eta(1:end-1,:)   + MAT.Eta(2:end,:) )/2); % crystal settling speed
segSis([1 end],:) = 0;
segSis(:,[1 end]) = sds*segSis(:,[2 end-1]);

segFes      = -2/9 .* ((MAT.rhoFes(1:end-1,:)+ MAT.rhoFes(2:end,:))./2      ...
                     -(MAT.rhot(1:end-1,:)  + MAT.rhot(2:end,:))/2)...
                  .*  PHY.gz*PHY.d^2.* (phiFes(1:end-1,:)+phiFes(2:end,:))./2  ...
                  ./ ((MAT.Eta(1:end-1,:)   + MAT.Eta(2:end,:) )/2); % crystal settling speed
segFes([1 end],:) = 0;
segFes(:,[1 end]) = sds*segFes(:,[2 end-1]);

segFel      = -2/9 .* ((MAT.rhoFel(1:end-1,:)+ MAT.rhoFel(2:end,:))./2 ...
                     -(MAT.rhot(1:end-1,:)  + MAT.rhot(2:end,:))/2)...
                  .*  PHY.gz*PHY.d^2.* (phiFel(1:end-1,:)+phiFel(2:end,:))./2 ...
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

% end

% update velocity divergence
Div_V(2:end-1,2:end-1) =...
    ddz(SOL.W(:,2:end-1),NUM.h) ...                           % get velocity divergence
                       + ddx(SOL.U(2:end-1,:),NUM.h);
Div_V([1 end],:) = Div_V([2 end-1],:);                                     % apply boundary conditions
Div_V(:,[1 end]) = Div_V(:,[2 end-1]);

% update volume source
Div_rhov =  + advection(MAT.rhoSis.*phiSis,0.*SOL.U,segSis,NUM.h,NUM.h,ADVN,'flx') ...
            + advection(MAT.rhoFes.*phiFes,0.*SOL.U,segFes,NUM.h,NUM.h,ADVN,'flx') ...
            + advection(MAT.rhoFel.*phiFel,0.*SOL.U,segFel,NUM.h,NUM.h,ADVN,'flx') ...
            + advection(MAT.rhot  ,   SOL.U,SOL.W ,NUM.h,NUM.h,ADVN,'adv');

VolSrc = -((MAT.rhot-MAT.rhoo)./NUM.dt + (Div_rhov - MAT.rhot.*Div_V + Div_rhoVo)/2)./(MAT.rhot/2);

dVoldt = mean(mean(VolSrc(2:end-1,2:end-1)));
VolSrc = VolSrc - dVoldt;

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