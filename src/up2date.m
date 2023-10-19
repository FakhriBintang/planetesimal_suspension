% print update header
% fprintf(1,'  ---  update materials & deformation \n');
tic;  % start clock on update

%% update T and C-dependent density
rhosSi = rhosSi0.*(1 - aT.*(T-perTSi) - gCSi.*(csSi-cphsSi1));
rholSi = rholSi0.*(1 - aT.*(T-perTSi) - gCSi.*(clSi-cphsSi1));
rhosFe = rhosFe0.*(1 - aT.*(T-perTFe) - gCFe.*(csFe-cphsFe1));
rholFe = rholFe0.*(1 - aT.*(T-perTFe) - gCFe.*(clFe-cphsFe1));

% update mixture density
rho    = 1./(xFe.*(fsFe./rhosFe + flFe./rholFe) ...
           + xSi.*(fsSi./rhosSi + flSi./rholSi));
rho([1 end],:) = rho([2 end-1],:);  rho(:,[1 end]) = rho(:,[2 end-1]);


%% convert weight to volume fractions
phisFe = max(0,min(1,xFe.* fsFe .* rho ./ rhosFe));
philFe = max(0,min(1,xFe.* flFe .* rho ./ rholFe));
phisSi = max(0,min(1,xSi.* fsSi .* rho ./ rhosSi));
philSi = max(0,min(1,xSi.* flSi .* rho ./ rholSi)); 


%% update viscosity
% get pure phase viscosity
% Simple Arrhenian model for a melt of fixed composition
etas   = zeros(size(phisSi)) + EtaSol0; 
etalSi = EtalSi0 .* exp(Em./(8.3145.*(T+273.15))-Em./(8.3145.*(perTSi+273.15)));
etalFe = EtalFe0 .* exp(Em./(8.3145.*(T+273.15))-Em./(8.3145.*(perTFe+273.15)));

% get permission weights
kv = permute(cat(3,etas,etalSi,etalFe),[3,1,2]);
Mv = permute(repmat(kv,1,1,1,3),[4,1,2,3])./permute(repmat(kv,1,1,1,3),[1,4,2,3]);
 
ff = max(1e-6,min(1-1e-6,permute(cat(3,phisSi+phisFe,philSi,philFe),[3,1,2])));
FF = permute(repmat(ff,1,1,1,3),[4,1,2,3]);
Sf = (FF./BBP).^(1./CCP);  Sf = Sf./sum(Sf,2);
Xf = sum(AAP.*Sf,2).*FF + (1-sum(AAP.*Sf,2)).*Sf;

% get momentum and volume permissions
thtv = squeeze(prod(Mv.^Xf,2));

% get momentum flux and transfer coefficients
Cv = ((1-ff)./[dx;dm;df].^2.*kv.*thtv).^-1;

% compose effective viscosity, segregation coefficients
Eta  = squeeze(sum(ff.*kv.*thtv,1));                                   % effective magma viscosity
etamax   = 1e+6.*min(Eta(:));                                          % set max eta for 1e6 max range
Eta  = (1./etamax + 1./(Eta * etareg)).^(-1);                      % limit viscosity range
Eta([1 end],:) = Eta([2 end-1],:);  
Eta(:,[1 end]) = Eta(:,[2 end-1]);
EtaC = (Eta(1:end-1,1:end-1).*Eta(2:end,1:end-1) ...           % viscosity in cell corners
         .* Eta(1:end-1,2:end  ).*Eta(2:end,2:end  )).^0.25;


%% update segregation cofficients
Ksgr_x = squeeze(Cv(1,:,:)) + TINY^2;  Ksgr_x([1 end],:) = Ksgr_x([2 end-1],:);  Ksgr_x(:,[1 end]) = Ksgr_x(:,[2 end-1]);
Ksgr_m = squeeze(Cv(2,:,:)) + TINY^2;  Ksgr_m([1 end],:) = Ksgr_m([2 end-1],:);  Ksgr_m(:,[1 end]) = Ksgr_m(:,[2 end-1]);
Ksgr_f = squeeze(Cv(3,:,:)) + TINY^2;  Ksgr_f([1 end],:) = Ksgr_f([2 end-1],:);  Ksgr_f(:,[1 end]) = Ksgr_f(:,[2 end-1]);
% dampen liquid segregation at low crystalinities
Ksgr_m = Ksgr_m.*(1-philSi).^2; Ksgr_m = max(TINY^2,Ksgr_m);
Ksgr_f = Ksgr_f.*(1-philFe).^2; Ksgr_f = max(TINY^2,Ksgr_f);


%% calculate stokes settling velocity
sds = BCsides;      % side boundary type 

% test alternative interpolations for segregation coefficients
segsSi = ((rhosSi(1:end-1,:)+rhosSi(2:end,:))./2-(rho(1:end-1,:)+rho(2:end,:))./2).*gz.*(Ksgr_x(1:end-1,:).*Ksgr_x(2:end,:)).^0.5;
segsFe = ((rhosFe(1:end-1,:)+rhosFe(2:end,:))./2-(rho(1:end-1,:)+rho(2:end,:))./2).*gz.*(Ksgr_x(1:end-1,:).*Ksgr_x(2:end,:)).^0.5;
seglSi = ((rholSi(1:end-1,:)+rholSi(2:end,:))./2-(rho(1:end-1,:)+rho(2:end,:))./2).*gz.*(Ksgr_m(1:end-1,:).*Ksgr_m(2:end,:)).^0.5; 
seglFe = ((rholFe(1:end-1,:)+rholFe(2:end,:))./2-(rho(1:end-1,:)+rho(2:end,:))./2).*gz.*(Ksgr_f(1:end-1,:).*Ksgr_f(2:end,:)).^0.5;

% zero boundary condition
segsSi([1 end],:) = 0;
segsSi(:,[1 end]) = sds*segsSi(:,[2 end-1]);
segsFe([1 end],:) = 0;
segsFe(:,[1 end]) = sds*segsFe(:,[2 end-1]);
seglSi([1 end],:) = 0;
seglSi(:,[1 end]) = sds*seglSi(:,[2 end-1]);
seglFe([1 end],:) = 0;
seglFe(:,[1 end]) = sds*seglFe(:,[2 end-1]);


%% update diffusion parameters
kT    =  (xFe.*kTFe + xSi.*kTSi);               % magma thermal conductivity
ks    = kT./T;
Vel(2:end-1,2:end-1) = sqrt(((W(1:end-1,2:end-1)+W(2:end,2:end-1))/2).^2 ...
                  + ((U(2:end-1,1:end-1)+U(2:end-1,2:end))/2).^2);
Vel([1 end],:) = Vel([2 end-1],:);
Vel(:,[1 end]) = Vel(:,[2 end-1]);
kW  = Vel/10*h/10; % diffusivity due to turbulent eddies 
kwlFe = abs((rholFe-rho).*gz0.*Ksgr_f*df*10);                                   % segregation fluctuation diffusivity
kwsFe = abs((rhosFe-rho).*gz0.*Ksgr_x*dx*10); 
kwsSi = abs((rhosSi-rho).*gz0.*Ksgr_x*dx*10); 
klFe  = philSi.*philFe.*(kwlFe + kW);
ksFe  = philSi.*phisFe.*(kwsFe + kW);
ksSi  = philSi.*phisSi.*(kwsSi + kW);



%% velocity divergence and volume source
% update velocity divergence
Div_V(2:end-1,2:end-1) = ddz(W(:,2:end-1),h) ...                   % get velocity divergence
                       + ddx(U(2:end-1,:),h);
Div_V([1 end],:) = Div_V([2 end-1],:);                                     % apply boundary conditions
Div_V(:,[1 end]) = Div_V(:,[2 end-1]);

% update volume source
if step>0
    Div_rhoV = + advect(FsSi(inz,inx),UsSi(inz,:),WsSi(:,inx),h,{ADVN,''   },[1,2],BCA) ...
               + advect(FlSi(inz,inx),UlSi(inz,:),WlSi(:,inx),h,{ADVN,''   },[1,2],BCA) ...
               + advect(FsFe(inz,inx),UsFe(inz,:),WsFe(:,inx),h,{ADVN,''   },[1,2],BCA) ...
               + advect(FlFe(inz,inx),UlFe(inz,:),WlFe(:,inx),h,{ADVN,''   },[1,2],BCA);
    F_DivV   = (alpha1*rho(inz,inx) - alpha2*rhoo(inz,inx) - alpha3*rhooo(inz,inx))./dt + (beta1*Div_rhoV + beta2*Div_rhoVo + beta3*Div_rhoVoo);  % get residual of mixture mass conservation
    VolSrc   = Div_V(inz,inx) - F_DivV./rho(inz,inx)/2;  % correct volume source term by scaled residual
end

% set variable boundary conditions
UBG    = mean(mean(VolSrc))./2 .* (L/2 - XU);
WBG    = mean(mean(VolSrc))./2 .* (D/2 - ZW);


%% Calculate stress and strain rates
% update strain rates
SOLexx(:,2:end-1) = diff(U,1,2)./h - Div_V(:,2:end-1)./3;                     % x-normal strain rate
SOLexx([1 end],:) = SOLexx([2 end-1],:);                                         % apply boundary conditions
SOLexx(:,[1 end]) = SOLexx(:,[2 end-1]);
SOLezz(2:end-1,:) = diff(W,1,1)./h - Div_V(2:end-1,:)./3;                        % z-normal strain rate
SOLezz([1 end],:) = SOLezz([2 end-1],:);                                         % apply boundary conditions
SOLezz(:,[1 end]) = SOLezz(:,[2 end-1]);
SOLexz            = 1/2.*(diff(U,1,1)./h+diff(W,1,2)./h);                     % shear strain rate

% update stresses
SOLtxx = Eta .* SOLexx;                                                          % x-normal stress
SOLtzz = Eta .* SOLezz;                                                          % z-normal stress
SOLtxz = EtaC.* SOLexz;   

% entropy production/heat dissipation rate
[grdTx,grdTz] = gradient(T,h);
EntProd = ks(2:end-1,2:end-1).*(grdTz(2:end-1,2:end-1).^2 + grdTx(2:end-1,2:end-1).^2)...
    + Hr(2:end-1,2:end-1) ... % Hr in W/m^3, if we instead choose W/Kg, miltiply by \rho
    + SOLexx(2:end-1,2:end-1).*SOLtxx(2:end-1,2:end-1) +SOLexx(2:end-1,2:end-1).*SOLtxx(2:end-1,2:end-1)...
    + 2.*(SOLexz(1:end-1,1:end-1)+SOLexz(2:end,1:end-1)+SOLexz(1:end-1,2:end)+SOLexz(2:end,2:end))./4 ...
    .*(SOLtxz(1:end-1,1:end-1)+SOLtxz(2:end,1:end-1)+SOLtxz(1:end-1,2:end)+SOLtxz(2:end,2:end))./4 ...
    +  phisFe(2:end-1,2:end-1)./Ksgr_x(2:end-1,2:end-1) .* ((segsFe(1:end-1,2:end-1)+segsFe(2:end,2:end-1))./2).^2  ...
    +  philFe(2:end-1,2:end-1)./Ksgr_f(2:end-1,2:end-1) .* ((seglFe(1:end-1,2:end-1)+seglFe(2:end,2:end-1))./2).^2  ...
    +  phisSi(2:end-1,2:end-1)./Ksgr_x(2:end-1,2:end-1) .* ((segsSi(1:end-1,2:end-1)+segsSi(2:end,2:end-1))./2).^2  ...
    +  philSi(2:end-1,2:end-1)./Ksgr_m(2:end-1,2:end-1) .* ((seglSi(1:end-1,2:end-1)+seglSi(2:end,2:end-1))./2).^2;

