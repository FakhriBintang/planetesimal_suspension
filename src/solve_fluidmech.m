Wi = SOL.W; Ui = SOL.U; Pi = SOL.P;
% get mapping arrays
NW = NUM.NW; NU = NUM.NU; NP = NUM.NP;
% profile on

%% assemble coefficients for matrix velocity diagonal and right-hand side
II  = [];       % equation indeces into A
JJ  = [];       % variable indeces into A
AA  = [];       % coefficients for A
IR  = [];       % equation indeces into R
RR  = [];       % forcing entries for R

% assemble coefficients of z-stress divergence
    
% left boundary
ii = NUM.MapW(:,1); jj1 = ii; jj2 = NUM.MapW(:,2);
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)+SOL.BCsides];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% right boundary
ii = NUM.MapW(:,end); jj1 = ii; jj2 = NUM.MapW(:,end-1);
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)+SOL.BCsides];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% top boundary
ii = NUM.MapW(1,2:end-1); jj = ii;
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj(:)];   AA = [AA; aa(:)+1];
aa = zeros(size(ii)) - WBG(1,2:end-1); 
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% bottom boundary
ii = NUM.MapW(end,2:end-1); jj = ii;
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj(:)];   AA = [AA; aa(:)+1];
aa = zeros(size(ii)) - WBG(end,2:end-1);
IR = [IR; ii(:)]; RR = [RR; aa(:)];


% internal points
ii    = NUM.MapW(2:end-1,2:end-1);
EtaC1 = MAT.EtaC(2:end-1,1:end-1); EtaC2 = MAT.EtaC(2:end-1,2:end  );
EtaP1 = MAT.Eta (2:end-2,2:end-1); EtaP2 = MAT.Eta (3:end-1,2:end-1);

% coefficients multiplying z-velocities W
%             top          ||         bottom          ||           left            ||          right
jj1 = NUM.MapW(1:end-2,2:end-1); jj2 = NUM.MapW(3:end,2:end-1); jj3 = NUM.MapW(2:end-1,1:end-2); jj4 = NUM.MapW(2:end-1,3:end);

aa = - 2/3*(EtaP1+EtaP2)/NUM.h^2 - 1/2*(EtaC1+EtaC2)/NUM.h^2;
II = [II; ii(:)]; JJ = [JJ;  ii(:)];   AA = [AA; aa(:)               ];      % W on stencil centre
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; 2/3*EtaP1(:)/NUM.h^2];      % W one above
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; 2/3*EtaP2(:)/NUM.h^2];      % W one below
II = [II; ii(:)]; JJ = [JJ; jj3(:)];   AA = [AA; 1/2*EtaC1(:)/NUM.h^2];      % W one to the left
II = [II; ii(:)]; JJ = [JJ; jj4(:)];   AA = [AA; 1/2*EtaC2(:)/NUM.h^2];      % W one to the right

% what shall we do with a drunken sailor...
aa = -ddz(MAT.rho(2:end-1,2:end-1),NUM.h).*MAT.gz(2:end-1,2:end-1).*NUM.dt/2;
II = [II; ii(:)]; JJ = [JJ;  ii(:)];   AA = [AA; aa(:)];

% coefficients multiplying x-velocities U
%         top left         ||        bottom left          ||       top right       ||       bottom right
jj1 = NUM.MapU(2:end-2,1:end-1); jj2 = NUM.MapU(3:end-1,1:end-1); jj3 = NUM.MapU(2:end-2,2:end); jj4 = NUM.MapU(3:end-1,2:end);

II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; (1/2*EtaC1(:)-1/3*EtaP1(:))/NUM.h/NUM.h];  % W one to the top and left
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA;-(1/2*EtaC1(:)-1/3*EtaP2(:))/NUM.h/NUM.h];  % W one to the bottom and left
II = [II; ii(:)]; JJ = [JJ; jj3(:)];   AA = [AA;-(1/2*EtaC2(:)-1/3*EtaP1(:))/NUM.h/NUM.h];  % W one to the top and right
II = [II; ii(:)]; JJ = [JJ; jj4(:)];   AA = [AA; (1/2*EtaC2(:)-1/3*EtaP2(:))/NUM.h/NUM.h];  % W one to the bottom and right


% z-RHS vector
if ~RUN.bnchm
    rhoBF =    (MAT.rho(2:end-2,2:end-1) + MAT.rho(3:end-1,2:end-1))/2 - rhoRef;
    if NUM.nxP<=10; rhoBF = repmat(mean(rhoBF,2),1,NUM.nxW-2); end
end

rr = - rhoBF .* MAT.gz(2:end-1,2:end-1);
if RUN.bnchm; rr = rr + src_W_mms(2:end-1,2:end-1); end


IR = [IR; ii(:)];  RR = [RR; rr(:)];


%  assemble coefficients of x-stress divergence

% top boundary
ii = NUM.MapU(1,:); jj1 = ii; jj2 = NUM.MapU(2,:);
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)+SOL.BCtop];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% bottom boundary
ii = NUM.MapU(end,:); jj1 = ii; jj2 = NUM.MapU(end-1,:);
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)+SOL.BCbot];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% left side boundary
ii = NUM.MapU(2:end-1,1); jj = ii;
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj(:)];   AA = [AA; aa(:)+1];
aa = zeros(size(ii)) - UBG(2:end-1,1);
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% right side boundary
ii = NUM.MapU(2:end-1,end); jj = ii;
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj(:)];   AA = [AA; aa(:)+1];
aa = zeros(size(ii)) - UBG(2:end-1,end);
IR = [IR; ii(:)]; RR = [RR; aa(:)];


% internal points
ii    = NUM.MapU(2:end-1,2:end-1);
EtaC1 = MAT.EtaC(1:end-1,2:end-1); EtaC2 = MAT.EtaC(2:end  ,2:end-1);
EtaP1 = MAT.Eta (2:end-1,2:end-2); EtaP2 = MAT.Eta (2:end-1,3:end-1);

% coefficients multiplying x-velocities U
%            left               ||          right             ||           top                 ||          bottom
jj1 = NUM.MapU(2:end-1,1:end-2); jj2 = NUM.MapU(2:end-1,3:end); jj3 = NUM.MapU(1:end-2,2:end-1); jj4 = NUM.MapU(3:end,2:end-1);

aa = - 2/3*(EtaP1+EtaP2)/NUM.h^2 - 1/2*(EtaC1+EtaC2)/NUM.h^2;
II = [II; ii(:)]; JJ = [JJ;  ii(:)];   AA = [AA; aa(:)           ];      % U on stencil centre
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; 2/3*EtaP1(:)/NUM.h^2];      % U one to the left
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; 2/3*EtaP2(:)/NUM.h^2];      % U one to the right
II = [II; ii(:)]; JJ = [JJ; jj3(:)];   AA = [AA; 1/2*EtaC1(:)/NUM.h^2];      % U one above
II = [II; ii(:)]; JJ = [JJ; jj4(:)];   AA = [AA; 1/2*EtaC2(:)/NUM.h^2];      % U one below

% coefficients multiplying z-velocities W
%         top left         ||        top right          ||       bottom left       ||       bottom right
jj1 = NUM.MapW(1:end-1,2:end-2); jj2 = NUM.MapW(1:end-1,3:end-1); jj3 = NUM.MapW(2:end,2:end-2); jj4 = NUM.MapW(2:end,3:end-1);

II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; (1/2*EtaC1(:)-1/3*EtaP1(:))/NUM.h/NUM.h];  % W one to the top and left
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA;-(1/2*EtaC1(:)-1/3*EtaP2(:))/NUM.h/NUM.h];  % W one to the top and right
II = [II; ii(:)]; JJ = [JJ; jj3(:)];   AA = [AA;-(1/2*EtaC2(:)-1/3*EtaP1(:))/NUM.h/NUM.h];  % W one to the bottom and left
II = [II; ii(:)]; JJ = [JJ; jj4(:)];   AA = [AA; (1/2*EtaC2(:)-1/3*EtaP2(:))/NUM.h/NUM.h];  % W one to the bottom and right


% x-RHS vector
rr = zeros(size(ii)); % no x-buoyancy
if RUN.bnchm; rr = rr + src_U_mms(2:end-1,2:end-1); end

IR = [IR; ii(:)];  RR = [RR; rr(:)];


% assemble coefficient matrix & right-hand side vector
KV = sparse(II,JJ,AA,NUM.NW+NUM.NU,NUM.NW+NUM.NU);
RV = sparse(IR,ones(size(IR)),RR);


%% assemble coefficients for gradient operator
% reset global lists for vectorised assembly
II  = [];   JJ  = [];   AA  = [];

% Z-Stokes equation
ii  = NUM.MapW(2:end-1,2:end-1);
%             top              ||          bottom
jj1 = NUM.MapP(2:end-2,2:end-1); jj2 = NUM.MapP(3:end-1,2:end-1);

aa  = zeros(size(ii));
II  = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)-1/NUM.h];     % P one to the top
II  = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)+1/NUM.h];     % P one to the bottom

% X-Stokes equation
ii  = NUM.MapU(2:end-1,2:end-1);
%             left             ||           right
jj1 = NUM.MapP(2:end-1,2:end-2); jj2 = NUM.MapP(2:end-1,3:end-1);

aa  = zeros(size(ii));
II  = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)-1/NUM.h];     % P one to the left
II  = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)+1/NUM.h];     % P one to the right

GG  = sparse(II,JJ,AA,NW+NU,NP);


%% assemble coefficients for divergence operator
II  = [];       % equation indeces into A
JJ  = [];       % variable indeces into A
AA  = [];       % coefficients for A

%internal points
ii = NUM.MapP(2:end-1,2:end-1);

% coefficients multiplying velocities U, W
%            left U            ||            right U          ||             top W             ||          bottom W
jj1 = NUM.MapU(2:end-1,1:end-1); jj2 = NUM.MapU(2:end-1,2:end); jj3 = NUM.MapW(1:end-1,2:end-1); jj4 = NUM.MapW(2:end,2:end-1);

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)-1/NUM.h];  % U one to the left
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)+1/NUM.h];  % U one to the right
II = [II; ii(:)]; JJ = [JJ; jj3(:)];   AA = [AA; aa(:)-1/NUM.h];  % W one above
II = [II; ii(:)]; JJ = [JJ; jj4(:)];   AA = [AA; aa(:)+1/NUM.h];  % W one below

% assemble coefficient matrix
DD = sparse(II,JJ,AA,NUM.NP,NUM.NW+NUM.NU);


%% assemble coefficients for matrix pressure diagonal and right-hand side
II  = [];       % equation indeces into A
JJ  = [];       % variable indeces into A
AA  = [];       % coefficients for A
IR  = [];       % equation indeces into R
RR  = [];       % forcing entries for R

% boundary points
ii  = [NUM.MapP(1,:).'; NUM.MapP(end  ,:).']; % top & bottom
jj1 = ii;
jj2 = [NUM.MapP(2,:).'; NUM.MapP(end-1,:).'];

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)-1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

ii  = [NUM.MapP(:,1); NUM.MapP(:,end  )]; % left & right
jj1 = ii;
jj2 = [NUM.MapP(:,2); NUM.MapP(:,end-1)];

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)-1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];


% internal points
ii = NUM.MapP(2:end-1,2:end-1);

% coefficients multiplying matrix pressure P
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; ii(:)];    AA = [AA; aa(:)];  % P on stencil centre

% RHS
rr = - VolSrc;
if RUN.bnchm; rr = rr + src_P_mms(2:end-1,2:end-1); end

IR = [IR; ii(:)];
RR = [RR; rr(:)];


% assemble coefficient matrix and right-hand side vector
KP = sparse(II,JJ,AA,NP,NP);
RP = sparse(IR,ones(size(IR)),RR,NP,1);

nzp = round((NUM.nzP-2)/2)+1;
nxp = round((NUM.nxP-2)/2)+1;
DD(NUM.MapP(nzp,nxp),:) = 0;
KP(NUM.MapP(nzp,nxp),:) = 0;
KP(NUM.MapP(nzp,nxp),NUM.MapP(nzp,nxp)) = 1;
RP(NUM.MapP(nzp,nxp),:) = 0;
if RUN.bnchm; RP(NUM.MapP(nzp,nxp),:) = P_mms(nzp,nxp); end

%% assemble global coefficient matrix and right-hand side vector
LL =  [ KV -GG ; ...
       -DD  KP ];

RR = [RV; RP];

SCL = sqrt(abs(diag(LL)));
SCL = diag(sparse(1./(SCL+1)));

LL  = SCL*LL*SCL;
RR  = SCL*RR;

%% Solve linear system of equations for vx, vz, P
S = SCL*(LL\RR);  % update solution

% Read out solution
% map solution vector to 2D arrays
SOL.W  = full(reshape(S(NUM.MapW(:)),NUM.nzW, NUM.nxW));         % matrix z-velocity
SOL.U  = full(reshape(S(NUM.MapU(:)),NUM.nzU, NUM.nxU));         % matrix x-velocity
SOL.P  = full(reshape(S(NUM.MapP(:)...
                    + NUM.NW+NUM.NU),NUM.nzP, NUM.nxP));         % matrix dynamic pressure

SOL.UP(:,2:end-1) = (SOL.U(:,1:end-1)+SOL.U(:,2:end))./2;
SOL.WP(2:end-1,:) = (SOL.W(1:end-1,:)+SOL.W(2:end,:))./2;

% get residual of fluid mechanics equations from iterative update
resnorm_VP = norm((SOL.W - Wi).*any(SOL.W(:)>1e-12),2)./(norm(SOL.W,2)+TINY) ...
           + norm((SOL.U - Ui).*any(SOL.U(:)>1e-12),2)./(norm(SOL.U,2)+TINY) ...
           + norm(SOL.P - Pi,2)./(norm(SOL.P,2)+TINY);


% update phase velocities
WlSi = SOL.W + MAT.seglSi;
UlSi = SOL.U;
WsSi = SOL.W + MAT.segsSi;
UsSi = SOL.U;
WlFe = SOL.W + MAT.seglFe;
UlFe = SOL.U;
WsFe = SOL.W + MAT.segsFe;
UsFe = SOL.U;


%% update physical time step
dtadvn =  NUM.h/2   /max(abs([UlSi(:);WlSi(:);UsSi(:);WsSi(:);UlFe(:);WlFe(:);UsFe(:);WsFe(:)])); % stable timestep for advection
dtdiff = (NUM.h/2)^2/max(max(MAT.ks(:).*SOL.T(:))./MAT.rho(:)./PHY.Cp);                         % stable time step for T diffusion

NUM.dt = min(min(dtdiff,NUM.CFL * dtadvn),dtmax);                      % fraction of minimum stable time step

if NUM.dt==dtmax
    dtlimit = 'max step limited';
elseif dtdiff<dtadvn
    dtlimit = 'diffusion limited';
else
    dtlimit = 'advection limited';
end

% profile report
% profile off
