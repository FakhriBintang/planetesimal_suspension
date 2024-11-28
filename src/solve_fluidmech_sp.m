% profile on

% update velocity
if step>0
res_rho   = (a1*rho(inz,inx) - a2*rhoo(inz,inx) - a3*rhooo(inz,inx))./dt - (b1*advn_RHO + b2*advn_RHOo + b3*advn_RHOoo);  % get residual of mixture mass conservation
% volume source and background velocity passed to fluid-mechanics solver
upd_rho = - alpha*res_rho./rho(2:end-1,2:end-1)/b1 + beta*upd_rho;
VolSrc  = Div_V(2:end-1,2:end-1) + upd_rho;  % correct volume source term by scaled residual
W(:,2) = -flipud([0; cumsum(flipud(VolSrc.*rp(inz,:).^2).*h)])./rw.^2;

% update pressure
Div_tz = ddz(rp.^2.*tzz(:,2:end-1),h)./rw.^2;           % z-stress divergence
rhoBF  =    (rho(1:end-1,2) + rho(2:end,2))/2 - rhoRef;
PrsSrc = - (rhoBF ) .* gz(:,2) - Div_tz;
P(:,2) = full(flipud(cumsum(flipud([PrsSrc;0])*h)));

W(:,1) = W(:,2); W(:,end) = W(:,2);
P(:,1) = P(:,2); P(:,end) = P(:,2);
end

UP(:,2:end-1) = (U(:,1:end-1)+U(:,2:end))./2;
WP(2:end-1,:) = (W(1:end-1,:)+W(2:end,:))./2;

resnorm_VP = 0;
if ~bnchm
    % set phase diffusion speeds
    grd_philSi_z = ddz(philSi,h);
    grd_philSi_x = ddx(philSi,h);

    grd_phisSi_z = ddz(phisSi,h);
    grd_phisSi_x = ddx(phisSi,h);

    grd_philFe_z = ddz(philFe,h);
    grd_philFe_x = ddx(philFe,h);

    grd_phisFe_z = ddz(phisFe,h);
    grd_phisFe_x = ddx(phisFe,h);

    klSiz = (klSi(1:end-1,:)+klSi(2:end,:))./2;
    klSix = (klSi(:,1:end-1)+klSi(:,2:end))./2;

    ksSiz = (ksSi(1:end-1,:)+ksSi(2:end,:))./2;
    ksSix = (ksSi(:,1:end-1)+ksSi(:,2:end))./2;

    klFez = (klFe(1:end-1,:)+klFe(2:end,:))./2;
    klFex = (klFe(:,1:end-1)+klFe(:,2:end))./2;

    ksFez = (ksFe(1:end-1,:)+ksFe(2:end,:))./2;
    ksFex = (ksFe(:,1:end-1)+ksFe(:,2:end))./2;

    sumkz = klSiz + ksSiz + klFez + ksFez;
    grd_phistar_z = klSiz./sumkz .* grd_philSi_z ...
        + ksSiz./sumkz .* grd_phisSi_z ...
        + klFez./sumkz .* grd_philFe_z ...
        + ksFez./sumkz .* grd_phisFe_z;

    sumkx = klSix + ksSix + klFex + ksFex;
    grd_phistar_x = klSix./sumkx .* grd_philSi_x ...
        + ksSix./sumkx .* grd_phisSi_x ...
        + klFex./sumkx .* grd_philFe_x ...
        + ksFex./sumkx .* grd_phisFe_x;

    qlSiz    = - klSiz .* (grd_philSi_z - grd_phistar_z);
    qlSix    = - klSix .* (grd_philSi_x - grd_phistar_x);

    qsSiz    = - ksSiz .* (grd_phisSi_z - grd_phistar_z);
    qsSix    = - ksSix .* (grd_phisSi_x - grd_phistar_x);

    qlFez    = - klFez .* (grd_philFe_z - grd_phistar_z);
    qlFex    = - klFex .* (grd_philFe_x - grd_phistar_x);

    qsFez    = - ksFez .* (grd_phisFe_z - grd_phistar_z);
    qsFex    = - ksFex .* (grd_phisFe_x - grd_phistar_x);

    wqlSi = qlSiz./max(1e-8,(philSi(1:end-1,:)+philSi(2:end,:))./2);
    uqlSi = qlSix./max(1e-8,(philSi(:,1:end-1)+philSi(:,2:end))./2);

    wqsSi = qsSiz./max(1e-8,(phisSi(1:end-1,:)+phisSi(2:end,:))./2);
    uqsSi = qsSix./max(1e-8,(phisSi(:,1:end-1)+phisSi(:,2:end))./2);

    wqlFe = qlFez./max(1e-8,(philFe(1:end-1,:)+philFe(2:end,:))./2);
    uqlFe = qlFex./max(1e-8,(philFe(:,1:end-1)+philFe(:,2:end))./2);

    wqsFe = qsFez./max(1e-8,(phisFe(1:end-1,:)+phisFe(2:end,:))./2);
    uqsFe = qsFex./max(1e-8,(phisFe(:,1:end-1)+phisFe(:,2:end))./2);

    % update phase velocities (reference + diffusion + segregation)
    WlSi = W + wqlSi + seglSi;
    UlSi = U + uqlSi;
    WsSi = W + wqsSi + segsSi;
    UsSi = U + uqsSi;
    WlFe = W + wqlFe + seglFe;
    UlFe = U + uqlFe;
    WsFe = W + wqsFe + segsFe;
    UsFe = U + uqsFe;


    %% update physical time step
    maxcmp = 0.01;
    dtadvn =  h/2   /max(abs([UlSi(:);WlSi(:);UsSi(:);WsSi(:);UlFe(:);WlFe(:);UsFe(:);WsFe(:)])); % stable timestep for advection
    % dtdiff = (h/2)^2/max(max(ks(:).*T(:))./rho(:)./Cp)*0.5;                         % stable time step for T diffusion
    dtdiff = (h/2)^2/max([kc(:);klSi(:);ksSi(:);klFe(:);ksFe(:);(kT(:)+ks(:).*T(:))./rho(:)./Cp(:)])/2 ;
    if step>0
    dtc = maxcmp./max(abs([advn_FsFe(:)./rho(2:end-1,2);advn_FlFe(:)./rho(2:end-1,2);advn_FsSi(:)./rho(2:end-1,2);advn_FlSi(:)./rho(2:end-1,2)]));
    else; dtc = dtadvn;
    end

    dt = min(1.01*dto,min(min([dtdiff,CFL * dtadvn,dtc]),dtmax));                      % fraction of minimum stable time step
    
    [~,mindtm] = min([dtdiff; CFL * dtadvn; dtc; dtmax]);

switch mindtm
    % case 1 ; dtlimit = 'step change limited';
    case 1 ; dtlimit = 'diffusion limited';
    case 2 ; dtlimit = 'advection limited';
    case 3 ; dtlimit = 'phase change limited';
    case 4 ; 'max step limited';
end
end
% profile report
% profile off
