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
II  = [];       % equation indeces into A
JJ  = [];       % variable indeces into A
AA  = [];       % coefficients for A
IR  = [];       % equation indeces into R
RR  = [];       % forcing entries for R

% assemble coefficients of z-stress divergence
    
% left boundary
ii = MapW(:,1); jj1 = ii; jj2 = MapW(:,2);
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)+BCsides];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% right boundary
ii = MapW(:,end); jj1 = ii; jj2 = MapW(:,end-1);
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)+BCsides];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% top boundary
ii = MapW(1,2:end-1); jj = ii;
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj(:)];   AA = [AA; aa(:)+1];
aa = zeros(size(ii)) - WBG(1,2:end-1); 
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% bottom boundary
ii = MapW(end,2:end-1); jj = ii;
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj(:)];   AA = [AA; aa(:)+1];
aa = zeros(size(ii)) - WBG(end,2:end-1);
IR = [IR; ii(:)]; RR = [RR; aa(:)];


% internal points
ii    = MapW(2:end-1,2:end-1);
EtaC1 = EtaC(2:end-1,1:end-1); EtaC2 = EtaC(2:end-1,2:end  );
EtaP1 = Eta (2:end-2,2:end-1); EtaP2 = Eta (3:end-1,2:end-1);

% coefficients multiplying z-velocities W
%             top          ||         bottom          ||           left            ||          right
jj1 = MapW(1:end-2,2:end-1); jj2 = MapW(3:end,2:end-1); jj3 = MapW(2:end-1,1:end-2); jj4 = MapW(2:end-1,3:end);

aa = - 2/3*(EtaP1+EtaP2)/h^2 - 1/2*(EtaC1+EtaC2)/h^2;
II = [II; ii(:)]; JJ = [JJ;  ii(:)];   AA = [AA; aa(:)               ];      % W on stencil centre
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; 2/3*EtaP1(:)/h^2];      % W one above
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; 2/3*EtaP2(:)/h^2];      % W one below
II = [II; ii(:)]; JJ = [JJ; jj3(:)];   AA = [AA; 1/2*EtaC1(:)/h^2];      % W one to the left
II = [II; ii(:)]; JJ = [JJ; jj4(:)];   AA = [AA; 1/2*EtaC2(:)/h^2];      % W one to the right

% what shall we do with a drunken sailor...
if ~bnchm
    aa = -ddz(rho(2:end-1,2:end-1),h).*gz(2:end-1,2:end-1).*dt/2;
    II = [II; ii(:)]; JJ = [JJ;  ii(:)];   AA = [AA; aa(:)];
end
% coefficients multiplying x-velocities U
%         top left         ||        bottom left          ||       top right       ||       bottom right
jj1 = MapU(2:end-2,1:end-1); jj2 = MapU(3:end-1,1:end-1); jj3 = MapU(2:end-2,2:end); jj4 = MapU(3:end-1,2:end);

II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; (1/2*EtaC1(:)-1/3*EtaP1(:))/h/h];  % W one to the top and left
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA;-(1/2*EtaC1(:)-1/3*EtaP2(:))/h/h];  % W one to the bottom and left
II = [II; ii(:)]; JJ = [JJ; jj3(:)];   AA = [AA;-(1/2*EtaC2(:)-1/3*EtaP1(:))/h/h];  % W one to the top and right
II = [II; ii(:)]; JJ = [JJ; jj4(:)];   AA = [AA; (1/2*EtaC2(:)-1/3*EtaP2(:))/h/h];  % W one to the bottom and right


% z-RHS vector
if ~bnchm
    rhoBF =    (rho(2:end-2,2:end-1) + rho(3:end-1,2:end-1))/2 - rhoRef;
    if nxP<=10; rhoBF = repmat(mean(rhoBF,2),1,nxW-2); end
end

% rr = - (rhoBF - mean(rhoBF,2)) .* gz(2:end-1,2:end-1);
rr = - (rhoBF ) .* gz(2:end-1,2:end-1);

if bnchm; rr = rr + src_W_mms(2:end-1,2:end-1); end


IR = [IR; ii(:)];  RR = [RR; rr(:)];


%  assemble coefficients of x-stress divergence

% top boundary
ii = MapU(1,:); jj1 = ii; jj2 = MapU(2,:);
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)+BCtop];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% bottom boundary
ii = MapU(end,:); jj1 = ii; jj2 = MapU(end-1,:);
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)+BCbot];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% left side boundary
ii = MapU(2:end-1,1); jj = ii;
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj(:)];   AA = [AA; aa(:)+1];
aa = zeros(size(ii)) - UBG(2:end-1,1);
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% right side boundary
ii = MapU(2:end-1,end); jj = ii;
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj(:)];   AA = [AA; aa(:)+1];
aa = zeros(size(ii)) - UBG(2:end-1,end);
IR = [IR; ii(:)]; RR = [RR; aa(:)];


% internal points
ii    = MapU(2:end-1,2:end-1);
EtaC1 = EtaC(1:end-1,2:end-1); EtaC2 = EtaC(2:end  ,2:end-1);
EtaP1 = Eta (2:end-1,2:end-2); EtaP2 = Eta (2:end-1,3:end-1);

% coefficients multiplying x-velocities U
%            left               ||          right             ||           top                 ||          bottom
jj1 = MapU(2:end-1,1:end-2); jj2 = MapU(2:end-1,3:end); jj3 = MapU(1:end-2,2:end-1); jj4 = MapU(3:end,2:end-1);

aa = - 2/3*(EtaP1+EtaP2)/h^2 - 1/2*(EtaC1+EtaC2)/h^2;
II = [II; ii(:)]; JJ = [JJ;  ii(:)];   AA = [AA; aa(:)           ];      % U on stencil centre
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; 2/3*EtaP1(:)/h^2];      % U one to the left
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; 2/3*EtaP2(:)/h^2];      % U one to the right
II = [II; ii(:)]; JJ = [JJ; jj3(:)];   AA = [AA; 1/2*EtaC1(:)/h^2];      % U one above
II = [II; ii(:)]; JJ = [JJ; jj4(:)];   AA = [AA; 1/2*EtaC2(:)/h^2];      % U one below

% coefficients multiplying z-velocities W
%         top left         ||        top right          ||       bottom left       ||       bottom right
jj1 = MapW(1:end-1,2:end-2); jj2 = MapW(1:end-1,3:end-1); jj3 = MapW(2:end,2:end-2); jj4 = MapW(2:end,3:end-1);

II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; (1/2*EtaC1(:)-1/3*EtaP1(:))/h/h];  % W one to the top and left
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA;-(1/2*EtaC1(:)-1/3*EtaP2(:))/h/h];  % W one to the top and right
II = [II; ii(:)]; JJ = [JJ; jj3(:)];   AA = [AA;-(1/2*EtaC2(:)-1/3*EtaP1(:))/h/h];  % W one to the bottom and left
II = [II; ii(:)]; JJ = [JJ; jj4(:)];   AA = [AA; (1/2*EtaC2(:)-1/3*EtaP2(:))/h/h];  % W one to the bottom and right


% x-RHS vector
rr = zeros(size(ii)); % no x-buoyancy
if bnchm; rr = rr + src_U_mms(2:end-1,2:end-1); end

IR = [IR; ii(:)];  RR = [RR; rr(:)];


% assemble coefficient matrix & right-hand side vector
KV = sparse(II,JJ,AA,NW+NU,NW+NU);
RV = sparse(IR,ones(size(IR)),RR);


%% assemble coefficients for gradient operator
% reset global lists for vectorised assembly
II  = [];   JJ  = [];   AA  = [];

% Z-Stokes equation
ii  = MapW(2:end-1,2:end-1);
%             top              ||          bottom
jj1 = MapP(2:end-2,2:end-1); jj2 = MapP(3:end-1,2:end-1);

aa  = zeros(size(ii));
II  = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)-1/h];     % P one to the top
II  = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)+1/h];     % P one to the bottom

% X-Stokes equation
ii  = MapU(2:end-1,2:end-1);
%             left             ||           right
jj1 = MapP(2:end-1,2:end-2); jj2 = MapP(2:end-1,3:end-1);

aa  = zeros(size(ii));
II  = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)-1/h];     % P one to the left
II  = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)+1/h];     % P one to the right

GG  = sparse(II,JJ,AA,NW+NU,NP);


%% assemble coefficients for divergence operator
II  = [];       % equation indeces into A
JJ  = [];       % variable indeces into A
AA  = [];       % coefficients for A

%internal points
ii = MapP(2:end-1,2:end-1);

% coefficients multiplying velocities U, W
%            left U            ||            right U          ||             top W             ||          bottom W
jj1 = MapU(2:end-1,1:end-1); jj2 = MapU(2:end-1,2:end); jj3 = MapW(1:end-1,2:end-1); jj4 = MapW(2:end,2:end-1);

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)-1/h];  % U one to the left
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)+1/h];  % U one to the right
II = [II; ii(:)]; JJ = [JJ; jj3(:)];   AA = [AA; aa(:)-1/h];  % W one above
II = [II; ii(:)]; JJ = [JJ; jj4(:)];   AA = [AA; aa(:)+1/h];  % W one below

% assemble coefficient matrix
DD = sparse(II,JJ,AA,NP,NW+NU);


%% assemble coefficients for matrix pressure diagonal and right-hand side
II  = [];       % equation indeces into A
JJ  = [];       % variable indeces into A
AA  = [];       % coefficients for A
IR  = [];       % equation indeces into R
RR  = [];       % forcing entries for R

% boundary points
ii  = [MapP(1,:).'; MapP(end  ,:).']; % top & bottom
jj1 = ii;
jj2 = [MapP(2,:).'; MapP(end-1,:).'];

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)-1];
% IR = [IR; ii(:)]; RR = [RR; aa(:)];

ii  = [MapP(:,1); MapP(:,end  )]; % left & right
jj1 = ii;
jj2 = [MapP(:,2); MapP(:,end-1)];

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)-1];
% IR = [IR; ii(:)]; RR = [RR; aa(:)];


% internal points
ii = MapP(2:end-1,2:end-1);

% coefficients multiplying matrix pressure P
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; ii(:)];    AA = [AA; aa(:)];  % P on stencil centre

% RHS
rr = - VolSrc;
if bnchm; rr = rr + src_P_mms(2:end-1,2:end-1); end

IR = [IR; ii(:)];
RR = [RR; rr(:)];


% assemble coefficient matrix and right-hand side vector
KP = sparse(II,JJ,AA,NP,NP);
RP = sparse(IR,ones(size(IR)),RR,NP,1);

nzp = round((nzP-2)/2)+1;
nxp = round((nxP-2)/2)+1;
DD(MapP(nzp,nxp),:) = 0;
KP(MapP(nzp,nxp),:) = 0;
KP(MapP(nzp,nxp),MapP(nzp,nxp)) = 1;
RP(MapP(nzp,nxp),:) = 0;
if bnchm; RP(MapP(nzp,nxp),:) = P_mms(nzp,nxp); end

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
W  = full(reshape(SOL(MapW(:)),nzW, nxW));         % matrix z-velocity
U  = full(reshape(SOL(MapU(:)),nzU, nxU));         % matrix x-velocity
P  = full(reshape(SOL(MapP(:)...
                    + NW+NU),nzP, nxP));         % matrix dynamic pressure

UP(:,2:end-1) = (U(:,1:end-1)+U(:,2:end))./2;
WP(2:end-1,:) = (W(1:end-1,:)+W(2:end,:))./2;

% get residual of fluid mechanics equations from iterative update
% resnorm_VP = norm((W - Wi).*any(W(:)>1e-12),2)./(norm(W,2)+TINY) ...
%            + norm((U - Ui).*any(U(:)>1e-12),2)./(norm(U,2)+TINY) ...
%            + norm(P - Pi,2)./(norm(P,2)+TINY);
resnorm_VP = 0;
if ~bnchm

    % set phase diffusion speeds
    qsSiz    = - (ksSi(1:end-1,:)+ksSi(2:end,:))./2 .* ddz(phisSi,h);
    qsSix    = - (ksSi(:,1:end-1)+ksSi(:,2:end))./2 .* ddx(phisSi,h);
    diff_sSi = - ddz(qsSiz(:,2:end-1),h)  ...                                 % heat diffusion
               - ddx(qsSix(2:end-1,:),h);
    qlFez    = - (klFe(1:end-1,:)+klFe(2:end,:))./2 .* ddz(philFe,h);
    qlFex    = - (klFe(:,1:end-1)+klFe(:,2:end))./2 .* ddx(philFe,h);
    diff_lFe = - ddz(qlFez(:,2:end-1),h)  ...                                 % heat diffusion
               - ddx(qlFex(2:end-1,:),h);
    qsFez    = - (ksFe(1:end-1,:)+ksFe(2:end,:))./2 .* ddz(phisFe,h);
    qsFex    = - (ksFe(:,1:end-1)+ksFe(:,2:end))./2 .* ddx(phisFe,h);
    diff_sFe = - ddz(qsFez(:,2:end-1),h)  ...                                 % heat diffusion
               - ddx(qsFex(2:end-1,:),h);

    qlSiz = qsSiz - qsFez - qlFez; qlSix = qsSix - qsFex - qlFex;
    
    wqlSi = qlSiz./max(1e-6,(philSi(1:end-1,:)+philSi(2:end,:))./2);
    uqlSi = qlSix./max(1e-6,(philSi(:,1:end-1)+philSi(:,2:end))./2);

    wqsSi = qsSiz./max(1e-6,(phisSi(1:end-1,:)+phisSi(2:end,:))./2);
    uqsSi = qsSix./max(1e-6,(phisSi(:,1:end-1)+phisSi(:,2:end))./2);

    wqlFe = qlFez./max(1e-6,(philFe(1:end-1,:)+philFe(2:end,:))./2);
    uqlFe = qlFex./max(1e-6,(philFe(:,1:end-1)+philFe(:,2:end))./2);

    wqsFe = qsFez./max(1e-6,(phisFe(1:end-1,:)+phisFe(2:end,:))./2);
    uqsFe = qsFex./max(1e-6,(phisFe(:,1:end-1)+phisFe(:,2:end))./2);
    
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
