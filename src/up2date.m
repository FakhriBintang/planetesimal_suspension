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

if selfgrav
[gz, gzP] = gravity(rho,rw,rho0);
end

%% convert weight to volume fractions
phisFe = max(0,min(1,xFe.* fsFe .* rho ./ rhosFe));
philFe = max(0,min(1,xFe.* flFe .* rho ./ rholFe));
phisSi = max(0,min(1,xSi.* fsSi .* rho ./ rhosSi));
philSi = max(0,min(1,xSi.* flSi .* rho ./ rholSi));


%% update viscosity
% get pure phase viscosity
% Simple Arrhenian model for a melt of fixed composition
etas   = ((etasSi.*phisSi)+(etasFe.*phisFe))./(phisSi+phisFe); % weighted average of solid viscosities
etas(isnan(etas)) = EtasSi0;
etalSi = EtalSi0 .* exp(Em./(8.3145.*T)-Em./(8.3145.*perTSi));
etalFe = EtalFe0 .* exp(Em./(8.3145.*T)-Em./(8.3145.*perTFe));

% update effective constituent sizes
dx = dx0.*(1-phisFe-phisSi).^dscale;
dm = dm0.*(1-philSi).^dscale;
df = df0.*(1-philFe).^dscale;

% get permission weights
kv = permute(cat(3,etas,etalSi,etalFe),[3,1,2]);
Mv = permute(repmat(kv,1,1,1,3),[4,1,2,3])./permute(repmat(kv,1,1,1,3),[1,4,2,3]);

dd = max(1e-6,min(1-1e-6,permute(cat(3,dx ,dm ,df ),[3,1,2])));
ff = max(1e-6,min(1-1e-6,permute(cat(3,phisSi+phisFe,philSi,philFe),[3,1,2])));
FF = permute(repmat(ff,1,1,1,3),[4,1,2,3]);
Sf = (FF./BBP).^(1./CCP);  Sf = Sf./sum(Sf,2);
Xf = sum(AAP.*Sf,2).*FF + (1-sum(AAP.*Sf,2)).*Sf;

% get momentum and volume permissions
thtv = squeeze(prod(Mv.^Xf,2));

% get momentum flux and transfer coefficients
Cv   = ff.*(1-ff)./dd.^2.*kv.*thtv;
Ksgr = ff./Cv;

% compose effective viscosity, segregation coefficients
Eta0  = squeeze(sum(ff.*kv.*thtv,1));                                   % effective magma viscosity


%% update segregation cofficients
Ksgr_x = squeeze(Ksgr(1,:,:)) + TINY^2;  Ksgr_x([1 end],:) = Ksgr_x([2 end-1],:);  Ksgr_x(:,[1 end]) = Ksgr_x(:,[2 end-1]); % crystal
Ksgr_m = squeeze(Ksgr(2,:,:)) + TINY^2;  Ksgr_m([1 end],:) = Ksgr_m([2 end-1],:);  Ksgr_m(:,[1 end]) = Ksgr_m(:,[2 end-1]); % silicate melt
Ksgr_f = squeeze(Ksgr(3,:,:)) + TINY^2;  Ksgr_f([1 end],:) = Ksgr_f([2 end-1],:);  Ksgr_f(:,[1 end]) = Ksgr_f(:,[2 end-1]); % iron melt

%% segregation coefficients for compaction length calculation
% Cv_m = squeeze(Cv(2,:,:)); % silicate melt
% Cv_f = squeeze(Cv(3,:,:)); % iron melt
% % dampened coefficients
% Cv_m = squeeze(ff(2,:,:))./Ksgr_m;
% Cv_f = squeeze(ff(3,:,:))./Ksgr_f;

%% calculate stokes settling velocity
sds = BCsides;      % side boundary type

% test alternative interpolations for segregation coefficients
segsSi = ((rhosSi(1:end-1,:)+rhosSi(2:end,:))./2-(rho(1:end-1,:)+rho(2:end,:))./2).*gz.*(Ksgr_x(1:end-1,:).*Ksgr_x(2:end,:)).^0.5;
segsFe = ((rhosFe(1:end-1,:)+rhosFe(2:end,:))./2-(rho(1:end-1,:)+rho(2:end,:))./2).*gz.*(Ksgr_x(1:end-1,:).*Ksgr_x(2:end,:)).^0.5;
seglSi = ((rholSi(1:end-1,:)+rholSi(2:end,:))./2-(rho(1:end-1,:)+rho(2:end,:))./2).*gz.*(Ksgr_m(1:end-1,:).*Ksgr_m(2:end,:)).^0.5;
seglFe = ((rholFe(1:end-1,:)+rholFe(2:end,:))./2-(rho(1:end-1,:)+rho(2:end,:))./2).*gz.*(Ksgr_f(1:end-1,:).*Ksgr_f(2:end,:)).^0.5;

if BCsegTop == 0
    % zero boundary (depletion) condition
    segsSi(1,:) = 0;
    segsFe(1,:) = 0;
    seglSi(1,:) = 0;
    seglFe(1,:) = 0;
else
    % continuous supply boundary condition
    segsSi(1,:) = segsSi(2,:);
    segsFe(1,:) = segsFe(2,:);
    seglSi(1,:) = seglSi(2,:);
    seglFe(1,:) = seglFe(2,:);
end

if BCsegBot == 0
    % zero boundary (accumulation) condition
    segsSi(end,:) = 0;
    segsFe(end,:) = 0;
    seglSi(end,:) = 0;
    seglFe(end,:) = 0;
else
    % sink boundary condition
    segsSi(end,:) = segsSi(end-1,:);
    segsFe(end,:) = segsFe(end-1,:);
    seglSi(end,:) = seglSi(end-1,:);
    seglFe(end,:) = seglFe(end-1,:);
end
segsSi(:,[1 end]) = sds*segsSi(:,[2 end-1]);
segsFe(:,[1 end]) = sds*segsFe(:,[2 end-1]);
seglSi(:,[1 end]) = sds*seglSi(:,[2 end-1]);
seglFe(:,[1 end]) = sds*seglFe(:,[2 end-1]);


%% velocity divergence and volume source
% update velocity divergence
switch mode
    case 'spherical'
    Div_V(2:end-1,2:end-1) = 1./(rp(2:end-1).^2).*ddz((rw.^2).*W(:,2),h);  % get velocity divergence
    Div_V([1 end],:) = Div_V([2 end-1],:);                                 % apply boundary conditions
    Div_V(:,[1 end]) = Div_V(:,[2 end-1]);
    otherwise
    Div_V(2:end-1,2:end-1) = ddz(W(:,2:end-1),h) ...                       % get velocity divergence
                           + ddx(U(2:end-1,:),h);
    Div_V([1 end],:) = Div_V([2 end-1],:);                                 % apply boundary conditions
    Div_V(:,[1 end]) = Div_V(:,[2 end-1]);
end


%% strain rates and velocity magnitude
% update strain rates
exx(:,2:end-1) = diff(U,1,2)./h - Div_V(:,2:end-1)./3;                     % x-normal strain rate
exx([1 end],:) = exx([2 end-1],:);                                         % apply boundary conditions
exx(:,[1 end]) = exx(:,[2 end-1]);
ezz(2:end-1,:) = diff(W,1,1)./h - Div_V(2:end-1,:)./3;                     % z-normal strain rate
ezz([1 end],:) = ezz([2 end-1],:);                                         % apply boundary conditions
ezz(:,[1 end]) = ezz(:,[2 end-1]);
exz            = 1/2.*(diff(U,1,1)./h+diff(W,1,2)./h);                     % shear strain rate

eII(2:end-1,2:end-1) = (0.5.*(exx(2:end-1,2:end-1).^2 + ezz(2:end-1,2:end-1).^2 ...
       + 2.*(exz(1:end-1,1:end-1).^2+exz(2:end,1:end-1).^2 ...
       +     exz(1:end-1,2:end  ).^2+exz(2:end,2:end  ).^2)/4)).^0.5 + TINY;
eII([1 end],:) = eII([2 end-1],:);
eII(:,[1 end]) = eII(:,[2 end-1]);

% velocity magnitude
Delta_cnv = h;
[~,drhodzz] = gradient(rho,h);
drhodz = max(0,-drhodzz); % + 1e-6.*rho;
Vel = drhodz.*gzP.*Delta_cnv.^3./Eta;

Vel([1 end],:) = Vel([2 end-1],:);
Vel(:,[1 end]) = Vel(:,[2 end-1]);


%% update diffusion parameters
if mixReg
    kW = Vel.*Delta_cnv;         % convective mixing diffusivity
else
    kW    = (h/2)^2.*eII + kmin;                                               % diffusivity due to turbulent eddies
end
    Pr = 3;
Sc = 3;

kT    = (xFe.*kTFe + xSi.*kTSi);                                           % magma thermal conductivity
% ks    = (kW+kT)./T;                                                             % entropy conductivity
ks    = (kW./Pr + kmin).*rho.*Cp./T;  
kc  =  kW./Sc + kmin;   
kwlFe = abs((rholFe-rho).*gzP.*Ksgr_f.*df*10) + kmin;                      % segregation fluctuation diffusivity
kwsFe = abs((rhosFe-rho).*gzP.*Ksgr_x.*dx*10) + kmin;
kwlSi = abs((rholSi-rho).*gzP.*Ksgr_m.*dm*10) + kmin;
kwsSi = abs((rhosSi-rho).*gzP.*Ksgr_x.*dx*10) + kmin;
klFe  = philFe.*(kwlFe+kW);
ksFe  = phisFe.*(kwsFe+kW);
klSi  = philSi.*(kwlSi+kW);
ksSi  = phisSi.*(kwsSi+kW);


%% Regularise and limit viscosity
Eta = (kW.*rho + Eta0)/2 + Eta/2;

etamax = 1e+8.*max(min(Eta(:)),etamin);                                                % set max eta for 1e8 max range
Eta    = 1./(1./etamax + 1./Eta) + etamin;                                       % limit viscosity range
% Eta([1 end],:) = Eta([2 end-1],:);
Eta(:,[1 end]) = Eta(:,[2 end-1]);
EtaC = (Eta(1:end-1,1:end-1).*Eta(2:end,1:end-1) ...           % viscosity in cell corners
     .* Eta(1:end-1,2:end  ).*Eta(2:end,2:end  )).^0.25;

%% Update dimensionless numbers
Ra     = Vel.*D/10./((kT+ks.*T) ./rho./Cp);
Re     = Vel.*D/10./( Eta       ./rho    );
%% update stresses
txx = Eta .* exx;                                                          % x-normal stress
tzz = Eta .* ezz;                                                          % z-normal stress
txz = EtaC.* exz;

tII(2:end-1,2:end-1) = (0.5.*(txx(2:end-1,2:end-1).^2 + tzz(2:end-1,2:end-1).^2 ...
       + 2.*(txz(1:end-1,1:end-1).^2+txz(2:end,1:end-1).^2 ...
       +     txz(1:end-1,2:end  ).^2+txz(2:end,2:end  ).^2)/4)).^0.5 + TINY;
tII([1 end],:) = tII([2 end-1],:);
tII(:,[1 end]) = tII(:,[2 end-1]);


%% entropy production/heat dissipation rate
[grdTx,grdTz] = gradient(T,h);
EntProd = ks(2:end-1,2:end-1).*(grdTz(2:end-1,2:end-1).^2 + grdTx(2:end-1,2:end-1).^2)...
    + exx(2:end-1,2:end-1).*txx(2:end-1,2:end-1) +exx(2:end-1,2:end-1).*txx(2:end-1,2:end-1)...
    + 2.*(exz(1:end-1,1:end-1)+exz(2:end,1:end-1)+exz(1:end-1,2:end)+exz(2:end,2:end))./4 ...
    .*(txz(1:end-1,1:end-1)+txz(2:end,1:end-1)+txz(1:end-1,2:end)+txz(2:end,2:end))./4 ...
    +  phisFe(2:end-1,2:end-1)./Ksgr_x(2:end-1,2:end-1) .* ((segsFe(1:end-1,2:end-1)+segsFe(2:end,2:end-1))./2).^2  ...
    +  philFe(2:end-1,2:end-1)./Ksgr_f(2:end-1,2:end-1) .* ((seglFe(1:end-1,2:end-1)+seglFe(2:end,2:end-1))./2).^2  ...
    +  phisSi(2:end-1,2:end-1)./Ksgr_x(2:end-1,2:end-1) .* ((segsSi(1:end-1,2:end-1)+segsSi(2:end,2:end-1))./2).^2  ...
    +  philSi(2:end-1,2:end-1)./Ksgr_m(2:end-1,2:end-1) .* ((seglSi(1:end-1,2:end-1)+seglSi(2:end,2:end-1))./2).^2;
