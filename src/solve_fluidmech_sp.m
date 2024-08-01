% Wi = W; Ui = U; Pi = P;
% profile on

if ~bnchm && step>0
    % update volume source
    res_rho   = (a1*rho(inz,inx) - a2*rhoo(inz,inx) - a3*rhooo(inz,inx))./dt - (b1*advn_RHO + b2*advn_RHOo + b3*advn_RHOoo);  % get residual of mixture mass conservation
    % volume source and background velocity passed to fluid-mechanics solver
    VolSrc  = Div_V(2:end-1,2:end-1) - alpha*res_rho./rho(2:end-1,2:end-1) + beta*upd_rho;  % correct volume source term by scaled residual
    upd_rho =       - alpha*res_rho./rho(2:end-1,2:end-1) + beta*upd_rho;
    % set variable boundary conditions
    UBG    = 0*mean(mean(VolSrc))./2 .* (L/2 - XU);
    WBG    = 2*mean(mean(VolSrc))./2 .* (D/2 - ZW);
end

%% assemble coefficients for matrix velocity diagonal and right-hand side
IIL  = [];       % equation indeces into A left
JJL  = [];       % variable indeces into A left
AAL  = [];       % coefficients for A left
IIR  = [];       % equation indeces into RHS
AAR  = [];       % forcing entries for RHS

% assemble coefficients of z-stress divergence
% top boundary
ii = MapW(1); jj = ii;
aa = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj(:)];   AAL = [AAL; aa(:)+1];
aa = zeros(size(ii)) - WBG(1,2); 
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% bottom boundary
ii = MapW(end); jj = ii;
aa = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj(:)];   AAL = [AAL; aa(:)+1];
aa = zeros(size(ii)) - WBG(end);
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];


% internal points
ii    = MapW(2:end-1);
%       above(-1)      ||        below(+1)/(0)
EP1 = Eta (2:end-2,2); EP2 = Eta (3:end-1,2);
RP1   = rp  (2:end-2);   RP2   = rp  (3:end-1);
% RW1   = rw  (1:end-2);   RW2   = rw  (3:end);
% internal
RWI   = rw (2:end-1);
 
% coefficients multiplying z-velocities W
%       top(-1)    ||     bottom(+1)(0)   ||%           left            ||          right
jj1 = MapW(1:end-2); jj2 = MapW(3:end); %jj3 = MapW(2:end-1,1:end-2); jj4 = MapW(2:end-1,3:end);

aa = - 2./(RWI.^2)./(h^2).*((RP2.^2).*EP2 + (RP1.^2).*EP1);
IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)               ];      % W on stencil centre
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; 2  ./(RWI.^2)./(h^2).*(RP1.^2).*EP1];     % W one above
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; 2  ./(RWI.^2)./(h^2).*(RP2.^2).*EP2];     % W one below
% II = [II; ii(:)]; JJ = [JJ; jj3(:)];   AAL = [AAL; 1/2*EtaC1(:)/h^2];      % W one to the left
% II = [II; ii(:)]; JJ = [JJ; jj4(:)];   AAL = [AAL; 1/2*EtaC2(:)/h^2];      % W one to the right

% % what shall we do with a drunken sailor...
% if ~bnchm
%     aa = -ddz(rho(2:end-1,2),h).*gz(2:end-1,2).*dt/2;
%     II = [II; ii(:)]; JJ = [JJ;  ii(:)];   AAL = [AAL; aa(:)];
% end
% z-RHS vector
if ~bnchm
    rhoBF =    (rho(2:end-2,2) + rho(3:end-1,2))/2 - rhoRef;
    if nxP<=10; rhoBF = repmat(mean(rhoBF,2),1,nxW-2); end
end

% rr = - (rhoBF - mean(rhoBF,2)) .* gz(2:end-1,2:end-1);
rr = - (rhoBF ) .* gz(2:end-1,2);

IIR = [IIR; ii(:)];  AAR = [AAR; rr(:)];

% assemble coefficient matrix & right-hand side vector
KV = sparse(IIL,JJL,AAL,NW,NW);
RV = sparse(IIR,ones(size(IIR)),AAR);


%% assemble coefficients for gradient operator
% reset global lists for vectorised assembly
IIL  = [];   JJL  = [];   AAL  = [];

% Z-Stokes equation
ii  = MapW(2:end-1);
%             top              ||          bottom
jj1 = MapP(2:end-2); jj2 = MapP(3:end-1);

%       above          ||       (below)
RW1   = rw  (1:end-1);   RW2   = rw  (2:end);
% internal
RPI   = rp (2:end-1);

aa  = zeros(size(ii));
IIL  = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)-1/h];     % P one to the top
IIL  = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+1/h];     % P one to the bottom

% % X-Stokes equation
% ii  = MapU(2:end-1,2:end-1);
% %             left             ||           right
% jj1 = MapP(2:end-1,2:end-2); jj2 = MapP(2:end-1,3:end-1);
% 
% aa  = zeros(size(ii));
% II  = [II; ii(:)]; JJ = [JJ; jj1(:)];   AAL = [AAL; aa(:)-1/h];     % P one to the left
% II  = [II; ii(:)]; JJ = [JJ; jj2(:)];   AAL = [AAL; aa(:)+1/h];     % P one to the right

GG  = sparse(IIL,JJL,AAL,NW,NP);


%% assemble coefficients for divergence operator
IIL  = [];       % equation indeces into A
JJL  = [];       % variable indeces into A
AAL  = [];       % coefficients for A

%internal points
ii = MapP(2:end-1);

% coefficients multiplying velocities U, W
%         top W             ||          bottom W
jj3 = MapW(1:end-1); jj4 = MapW(2:end);

aa = zeros(size(ii));
% II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AAL = [AAL; aa(:)-1/h];  % U one to the left
% II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AAL = [AAL; aa(:)+1/h];  % U one to the right
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL; aa(:)-1./(RPI.^2)/h.*(RW1.^2)];  % W one above
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL; aa(:)+1./(RPI.^2)/h.*(RW2.^2)];  % W one below

% assemble coefficient matrix
DD = sparse(IIL,JJL,AAL,NP,NW);


%% assemble coefficients for matrix pressure diagonal and right-hand side
IIL  = [];       % equation indeces into A
JJL  = [];       % variable indeces into A
AAL  = [];       % coefficients for A
IIR  = [];       % equation indeces into R
AAR  = [];       % forcing entries for R

% boundary points
ii  = [MapP(1).'; MapP(end  ).']; % top & bottom
jj1 = ii;
jj2 = [MapP(2).'; MapP(end-1).'];

aa = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];

% internal points
ii = MapP(2:end-1);

% coefficients multiplying matrix pressure P
aa = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; ii(:)];    AAL = [AAL; aa(:)];  % P on stencil centre

% RHS
rr = - VolSrc;
if bnchm; rr = rr + src_P_mms(2:end-1); end

IIR = [IIR; ii(:)];
AAR = [AAR; rr(:)];


% assemble coefficient matrix and right-hand side vector
KP = sparse(IIL,JJL,AAL,NP,NP);
RP = sparse(IIR,ones(size(IIR)),AAR,NP,1);

nzp = round((nzP-2)/2)+1;
nxp = round((nxP-2)/2)+1;
DD(MapP(nzp),:) = 0;
KP(MapP(nzp),:) = 0;
KP(MapP(nzp),MapP(nzp)) = 1;
RP(MapP(nzp),:) = 0;
if bnchm; RP(MapP(nzp),:) = P_mms(nzp); end

%% assemble global coefficient matrix and right-hand side vector
LL =  [ KV -GG ; ...
       -DD  KP ];

RR = [RV; RP];

SCL = sqrt(abs(diag(LL)));
SCL = diag(sparse(1./(SCL+1)));

LL  = SCL*LL*SCL;
RR  = SCL*RR;

%% Solve linear system of equations for vx, vz, P
SOL = SCL*(LL\RR);  % update solution


% Read out solution
% map solution vector to 2D arrays
W(:,2)  = full(reshape(SOL(MapW(:)),nzW,1));         % matrix z-velocity
% U  = full(reshape(SOL(MapU(:)),nzU, nxU));         % matrix x-velocity
P(:,2)  = full(reshape(SOL(MapP(:)+ NW),nzP,1));         % matrix dynamic pressure

W(:,1) = W(:,2); W(:,end) = W(:,2);
P(:,1) = P(:,2); P(:,end) = P(:,2);
UP(:,2:end-1) = (U(:,1:end-1)+U(:,2:end))./2;
WP(2:end-1,:) = (W(1:end-1,:)+W(2:end,:))./2;

% get residual of fluid mechanics equations from iterative update
% resnorm_VP = norm((W - Wi).*any(W(:)>1e-12),2)./(norm(W,2)+TINY) ...
%            + norm((U - Ui).*any(U(:)>1e-12),2)./(norm(U,2)+TINY) ...
%            + norm(P - Pi,2)./(norm(P,2)+TINY);
resnorm_VP = 0;
if ~bnchm
    if step>0
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

    else
        WlSi = W;
        UlSi = U;
        WsSi = W;
        UsSi = U;
        WlFe = W;
        UlFe = U;
        WsFe = W;
        UsFe = U;
    end

    %% update physical time step
    dtadvn =  h/2   /max(abs([UlSi(:);WlSi(:);UsSi(:);WsSi(:);UlFe(:);WlFe(:);UsFe(:);WsFe(:)])); % stable timestep for advection
    dtdiff = (h/2)^2/max(max(ks(:).*T(:))./rho(:)./Cp);                         % stable time step for T diffusion

    dt = min(min(dtdiff,CFL * dtadvn),dtmax);                      % fraction of minimum stable time step

    if dt==dtmax
        dtlimit = 'max step limited';
    elseif dtdiff<dtadvn
        dtlimit = 'diffusion limited';
    else
        dtlimit = 'advection limited';
    end
end
% profile report
% profile off
