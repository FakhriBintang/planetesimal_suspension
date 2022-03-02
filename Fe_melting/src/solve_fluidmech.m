% get mapping arrays
indU = NUM.MapU;
indW = NUM.MapW;
indP = NUM.MapP-NUM.NW-NUM.NU;

NW = NUM.NW; NU = NUM.NU; NP = NUM.NP;
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

% coefficients multiplying x-velocities U
%         top left         ||        bottom left          ||       top right       ||       bottom right
jj1 = NUM.MapU(2:end-2,1:end-1); jj2 = NUM.MapU(3:end-1,1:end-1); jj3 = NUM.MapU(2:end-2,2:end); jj4 = NUM.MapU(3:end-1,2:end);

II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; (1/2*EtaC1(:)-1/3*EtaP1(:))/NUM.h/NUM.h];  % W one to the top and left
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA;-(1/2*EtaC1(:)-1/3*EtaP2(:))/NUM.h/NUM.h];  % W one to the bottom and left
II = [II; ii(:)]; JJ = [JJ; jj3(:)];   AA = [AA;-(1/2*EtaC2(:)-1/3*EtaP1(:))/NUM.h/NUM.h];  % W one to the top and right
II = [II; ii(:)]; JJ = [JJ; jj4(:)];   AA = [AA; (1/2*EtaC2(:)-1/3*EtaP2(:))/NUM.h/NUM.h];  % W one to the bottom and right


% z-RHS vector
if ~RUN.bnchm
    rhoBF =    NUM.theta .*(MAT.rho(2:end-2,2:end-1)+MAT.rho(3:end-1,2:end-1))/2 ...
          + (1-NUM.theta).*(   rhoo(2:end-2,2:end-1)+   rhoo(3:end-1,2:end-1))/2 - rhoRef;
    if NUM.nxP<=10; rhoBF = repmat(mean(rhoBF,2),1,NUM.nxW-2); end
end

rr = - rhoBF .* PHY.gz(2:end-1,2:end-1);
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
ii  = indW(2:end-1,2:end-1);
%         top              ||          bottom
jj1 = indP(2:end-2,2:end-1); jj2 = indP(3:end-1,2:end-1);

aa  = zeros(size(ii));
II  = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)-1/NUM.h];     % P one to the top
II  = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)+1/NUM.h];     % P one to the bottom

% X-Stokes equation
ii  = indU(2:end-1,2:end-1);
%         left             ||           right
jj1 = indP(2:end-1,2:end-2); jj2 = indP(2:end-1,3:end-1);

aa  = zeros(size(ii));
II  = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)-1/NUM.h];     % P one to the left
II  = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)+1/NUM.h];     % P one to the right

GG  = sparse(II,JJ,AA,NW+NU,NP);


%% assemble coefficients for divergence operator
II  = [];       % equation indeces into A
JJ  = [];       % variable indeces into A
AA  = [];       % coefficients for A

%internal points
ii = indP(2:end-1,2:end-1);

% coefficients multiplying velocities U, W
%          left U          ||           right U       ||           top W           ||          bottom W
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
ii  = [indP(1,:).'; indP(end  ,:).']; % top & bottom
jj1 = ii;
jj2 = [indP(2,:).'; indP(end-1,:).'];

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)-1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

ii  = [indP(:,1); indP(:,end  )]; % left & right
jj1 = ii;
jj2 = [indP(:,2); indP(:,end-1)];

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)-1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];


% internal points
ii = indP(2:end-1,2:end-1);

% coefficients multiplying matrix pressure P
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; ii(:)];    AA = [AA; aa(:)];  % P on stencil centre

rr = - VolSrc(2:end-1,2:end-1);
if RUN.bnchm; rr = rr + src_P_mms(2:end-1,2:end-1); end

IR = [IR; ii(:)];
RR = [RR; rr(:)];


% assemble coefficient matrix and right-hand side vector
KP = sparse(II,JJ,AA,NUM.NP,NUM.NP);
RP = sparse(IR,ones(size(IR)),RR,NP,1);

Pscale = sqrt(geomean(MAT.Eta(:))/NUM.h^2);

nzp = round((NUM.nzP-2)/2)+1;
nxp = round((NUM.nxP-2)/2)+1;
KP(indP(nzp,nxp),:) = 0;
KP(indP(nzp,nxp),indP(nzp,nxp)) = Pscale;
RP(indP(nzp,nxp),:) = 0;
if RUN.bnchm; RP(indP(nzp,nxp),:) = P_mms(nzp,nxp); end

%% assemble global coefficient matrix and right-hand side vector
LL =  [ KV         -Pscale.*GG ; ...
       -Pscale.*DD  Pscale.*KP ];

RR = [RV; RP.*Pscale];


%% get residual
FF         = LL*S - RR;
resnorm_VP = norm(FF(:),2)./sqrt(length(S(:)));

% map residual vector to 2D arrays
res_W  = full(reshape(FF(NUM.MapW(:)),NUM.nzW,NUM.nxW));   % z-velocity residual
res_U  = full(reshape(FF(NUM.MapU(:)),NUM.nzU,NUM.nxU));   % x-velocity residual
res_P  = full(reshape(FF(NUM.MapP(:)),NUM.nzP,NUM.nxP));   % dynamic pressure residual


%% Solve linear system of equations for vx, vz, P
S = LL\RR;  % update solution

% Read out solution
% map solution vector to 2D arrays
SOL.W  = full(reshape(S(NUM.MapW(:)),NUM.nzW, NUM.nxW));         % matrix z-velocity
SOL.U  = full(reshape(S(NUM.MapU(:)),NUM.nzU, NUM.nxU));         % matrix x-velocity
SOL.P  = full(reshape(S(NUM.MapP(:)),NUM.nzP,NUM.nxP)).*Pscale;  % matrix dynamic pressure

SOL.UP(:,2:end-1) = (SOL.U(:,1:end-1)+SOL.U(:,2:end))./2;
SOL.WP(2:end-1,:) = (SOL.W(1:end-1,:)+SOL.W(2:end,:))./2;
