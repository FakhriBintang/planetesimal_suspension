% print update header
% fprintf(1,'  ---  update materials & deformation \n');
tic;  % start clock on update

%% update T and C-dependent density
MAT.rhosSi = PHY.rhosSi.*(1 - PHY.aT.*(SOL.T-CHM.perTSi) - PHY.gCSi.*(CHM.csSi-CHM.cphsSi1));
MAT.rholSi = PHY.rholSi.*(1 - PHY.aT.*(SOL.T-CHM.perTSi) - PHY.gCSi.*(CHM.clSi-CHM.cphsSi1));
MAT.rhosFe = PHY.rhosFe.*(1 - PHY.aT.*(SOL.T-CHM.perTFe) - PHY.gCFe.*(CHM.csFe-CHM.cphsFe1));
MAT.rholFe = PHY.rholFe.*(1 - PHY.aT.*(SOL.T-CHM.perTFe) - PHY.gCFe.*(CHM.clFe-CHM.cphsFe1));

% update mixture density
MAT.rho    = 1./(CHM.xFe.*CHM.fsFe./MAT.rhosFe + CHM.xFe.*CHM.flFe./MAT.rholFe ...
               + CHM.xSi.*CHM.fsSi./MAT.rhosSi + CHM.xSi.*CHM.flSi./MAT.rholSi);
MAT.rho([1 end],:) = MAT.rho([2 end-1],:);  MAT.rho(:,[1 end]) = MAT.rho(:,[2 end-1]);


%% convert weight to volume fractions
MAT.phisFe = max(0,min(1,CHM.xFe.* CHM.fsFe .* MAT.rho ./ MAT.rhosFe));
MAT.philFe = max(0,min(1,CHM.xFe.* CHM.flFe .* MAT.rho ./ MAT.rholFe));
MAT.phisSi = max(0,min(1,CHM.xSi.* CHM.fsSi .* MAT.rho ./ MAT.rhosSi));
MAT.philSi = max(0,min(1,CHM.xSi.* CHM.flSi .* MAT.rho ./ MAT.rholSi)); 


%% update viscosity
% get pure phase densities
etas   = zeros(size(MAT.phisSi)) + PHY.EtaSol0; 
etalSi = PHY.EtalSi0 .* exp(Em./(8.3145.*(SOL.T+273.15))-Em./(8.3145.*(CHM.perTSi+273.15)));
etalFe = PHY.EtalFe0 .* exp(Em./(8.3145.*(SOL.T+273.15))-Em./(8.3145.*(CHM.perTFe+273.15)));

% get permission weights
kv = permute(cat(3,etas,etalSi,etalFe),[3,1,2]);
Mv = permute(repmat(kv,1,1,1,3),[4,1,2,3])./permute(repmat(kv,1,1,1,3),[1,4,2,3]);
 
ff = max(1e-6,min(1-1e-6,permute(cat(3,MAT.phisSi+MAT.phisFe,MAT.philSi,MAT.philFe),[3,1,2])));
FF = permute(repmat(ff,1,1,1,3),[4,1,2,3]);
Sf = (FF./BBP).^(1./CCP);  Sf = Sf./sum(Sf,2);
Xf = sum(AAP.*Sf,2).*FF + (1-sum(AAP.*Sf,2)).*Sf;

% get momentum and volume permissions
thtv = squeeze(prod(Mv.^Xf,2));

% get momentum flux and transfer coefficients
Cv = ((1-ff)./[PHY.dx;PHY.dm;PHY.df].^2.*kv.*thtv).^-1;

% compose effective viscosity, segregation coefficients
MAT.Eta  = squeeze(sum(ff.*kv.*thtv,1));                                             % effective magma viscosity
MAT.Eta  = (1./etamax + 1./(MAT.Eta * etareg)).^(-1);                       % limit viscosity range
MAT.Eta([1 end],:) = MAT.Eta([2 end-1],:);  
MAT.Eta(:,[1 end]) = MAT.Eta(:,[2 end-1]);
MAT.EtaC = (MAT.Eta(1:end-1,1:end-1)+MAT.Eta(2:end,1:end-1) ...            % viscosity in cell corners
         +  MAT.Eta(1:end-1,2:end  )+MAT.Eta(2:end,2:end  ))./4;

%% update segregation cofficients
Ksgr_x = squeeze(Cv(1,:,:)) + TINY^2;  Ksgr_x([1 end],:) = Ksgr_x([2 end-1],:);  Ksgr_x(:,[1 end]) = Ksgr_x(:,[2 end-1]);
Ksgr_m = squeeze(Cv(2,:,:)) + TINY^2;  Ksgr_m([1 end],:) = Ksgr_m([2 end-1],:);  Ksgr_m(:,[1 end]) = Ksgr_m(:,[2 end-1]);
Ksgr_f = squeeze(Cv(3,:,:)) + TINY^2;  Ksgr_f([1 end],:) = Ksgr_f([2 end-1],:);  Ksgr_f(:,[1 end]) = Ksgr_f(:,[2 end-1]);
% dampen liquid segregation at low crystalinities
Ksgr_m = Ksgr_m.*(1-MAT.philSi).^2; Ksgr_m = max(TINY^2,Ksgr_m);
Ksgr_f = Ksgr_f.*(1-MAT.philFe).^2; Ksgr_f = max(TINY^2,Ksgr_f);


%% calculate stokes settling velocity
sds = SOL.BCsides;      % side boundary type 

% also runge kutta advection time stepping
kappaseg = 500;

% test alternative interpolations for segregation coefficients
% minimum  segregation coefficient
% segSis = ((MAT.rhoSis(1:end-1,:)+MAT.rhoSis(2:end,:))./2-(MAT.rho(1:end-1,:)+MAT.rho(2:end,:))./2).*PHY.gz.*min(Ksgr_x(1:end-1,:),Ksgr_x(2:end,:));
% segFes = ((MAT.rhoFes(1:end-1,:)+MAT.rhoFes(2:end,:))./2-(MAT.rho(1:end-1,:)+MAT.rho(2:end,:))./2).*PHY.gz.*min(Ksgr_x(1:end-1,:),Ksgr_x(2:end,:));
% 
% segSil = ((MAT.rhoSil(1:end-1,:)+MAT.rhoSil(2:end,:))./2-(MAT.rho(1:end-1,:)+MAT.rho(2:end,:))./2).*PHY.gz.*min(Ksgr_m(1:end-1,:),Ksgr_m(2:end,:)); 
% segFel = ((MAT.rhoFel(1:end-1,:)+MAT.rhoFel(2:end,:))./2-(MAT.rho(1:end-1,:)+MAT.rho(2:end,:))./2).*PHY.gz.*min(Ksgr_f(1:end-1,:),Ksgr_f(2:end,:)); 

% harmonic average
MAT.segsSi = ((MAT.rhosSi(1:end-1,:)+MAT.rhosSi(2:end,:))./2-(MAT.rho(1:end-1,:)+MAT.rho(2:end,:))./2).*PHY.gz.*2./(1./Ksgr_x(1:end-1,:)+1./Ksgr_x(2:end,:));
MAT.segsFe = ((MAT.rhosFe(1:end-1,:)+MAT.rhosFe(2:end,:))./2-(MAT.rho(1:end-1,:)+MAT.rho(2:end,:))./2).*PHY.gz.*2./(1./Ksgr_x(1:end-1,:)+1./Ksgr_x(2:end,:));
MAT.seglSi = ((MAT.rholSi(1:end-1,:)+MAT.rholSi(2:end,:))./2-(MAT.rho(1:end-1,:)+MAT.rho(2:end,:))./2).*PHY.gz.*2./(1./Ksgr_m(1:end-1,:)+1./Ksgr_m(2:end,:)); 
MAT.seglFe = ((MAT.rholFe(1:end-1,:)+MAT.rholFe(2:end,:))./2-(MAT.rho(1:end-1,:)+MAT.rho(2:end,:))./2).*PHY.gz.*2./(1./Ksgr_f(1:end-1,:)+1./Ksgr_f(2:end,:));

% zero boundary condition
MAT.segsSi([1 end],:) = 0;
MAT.segsSi(:,[1 end]) = sds*MAT.segsSi(:,[2 end-1]);
MAT.segsFe([1 end],:) = 0;
MAT.segsFe(:,[1 end]) = sds*MAT.segsFe(:,[2 end-1]);
MAT.seglSi([1 end],:) = 0;
MAT.seglSi(:,[1 end]) = sds*MAT.seglSi(:,[2 end-1]);
MAT.seglFe([1 end],:) = 0;
MAT.seglFe(:,[1 end]) = sds*MAT.seglFe(:,[2 end-1]);


%% update diffusion parameters
MAT.ks    =  (CHM.xFe.*PHY.kTFe + CHM.xSi.*PHY.kTSi)./SOL.T;               % magma thermal conductivity

% magma velocity magnitude
% Vel  = sqrt(((SOL.W([1,1:end],:)+SOL.W([1:end,end],:))/2).^2 ...
%           + ((SOL.U(:,[1,1:end])+SOL.U(:,[1:end,end]))/2).^2);
% kW   = Vel*NUM.h/100;
% kwFe = Vel.*abs((MAT.rhoFes-MAT.rho).*PHY.gz.*Ksgr_x.*PHY.dx);
% kwSi = Vel.*abs((MAT.rhoSis-MAT.rho).*PHY.gz.*Ksgr_x.*PHY.dx);

MAT.Ds = CHM.xFe.*CHM.fsFe.*CHM.dEntrFe + CHM.xSi.*CHM.fsSi.*CHM.dEntrSi;


%% velocity divergence and volume source
% update velocity divergence
Div_V(2:end-1,2:end-1) = ddz(SOL.W(:,2:end-1),NUM.h) ...                   % get velocity divergence
                       + ddx(SOL.U(2:end-1,:),NUM.h);
Div_V([1 end],:) = Div_V([2 end-1],:);                                     % apply boundary conditions
Div_V(:,[1 end]) = Div_V(:,[2 end-1]);

% update volume source
if NUM.step>0
    Div_rhoV = + advect(CHM.FsSi(inz,inx),UsSi(inz,:),WsSi(:,inx),NUM.h,{ADVN,''   },[1,2],BCA) ...
               + advect(CHM.FlSi(inz,inx),UlSi(inz,:),WlSi(:,inx),NUM.h,{ADVN,''   },[1,2],BCA) ...
               + advect(CHM.FsFe(inz,inx),UsFe(inz,:),WsFe(:,inx),NUM.h,{ADVN,''   },[1,2],BCA) ...
               + advect(CHM.FlFe(inz,inx),UlFe(inz,:),WlFe(:,inx),NUM.h,{ADVN,''   },[1,2],BCA);
    F_DivV   = (MAT.rho(inz,inx)-rhoo(inz,inx))./NUM.dt + NUM.theta.*Div_rhoV + (1-NUM.theta).*Div_rhoVo;  % get residual of mixture mass conservation
    VolSrc   = Div_V(inz,inx) - F_DivV./MAT.rho(inz,inx);  % correct volume source term by scaled residual
end

% set variable boundary conditions
UBG    = mean(mean(VolSrc))./2 .* (NUM.L/2 - NUM.XU);
WBG    = mean(mean(VolSrc))./2 .* (NUM.D/2 - NUM.ZW);


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
        + DEF.exx(2:end-1,2:end-1).*DEF.txx(2:end-1,2:end-1) +DEF.exx(2:end-1,2:end-1).*DEF.txx(2:end-1,2:end-1)...
        + 2.*(DEF.exz(1:end-1,1:end-1)+DEF.exz(2:end,1:end-1)+DEF.exz(1:end-1,2:end)+DEF.exz(2:end,2:end))./4 ...
           .*(DEF.txz(1:end-1,1:end-1)+DEF.txz(2:end,1:end-1)+DEF.txz(1:end-1,2:end)+DEF.txz(2:end,2:end))./4 ...
        +  MAT.phisFe(2:end-1,2:end-1)./Ksgr_x(2:end-1,2:end-1) .* ((MAT.segsFe(1:end-1,2:end-1)+MAT.segsFe(2:end,2:end-1))./2).^2  ...
        +  MAT.philFe(2:end-1,2:end-1)./Ksgr_f(2:end-1,2:end-1) .* ((MAT.seglFe(1:end-1,2:end-1)+MAT.seglFe(2:end,2:end-1))./2).^2  ...
        +  MAT.phisSi(2:end-1,2:end-1)./Ksgr_x(2:end-1,2:end-1) .* ((MAT.segsSi(1:end-1,2:end-1)+MAT.segsSi(2:end,2:end-1))./2).^2  ...
        +  MAT.philSi(2:end-1,2:end-1)./Ksgr_m(2:end-1,2:end-1) .* ((MAT.seglSi(1:end-1,2:end-1)+MAT.seglSi(2:end,2:end-1))./2).^2;

