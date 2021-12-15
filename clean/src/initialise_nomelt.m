% suspension flow: initialise model run

% print initialisation header
fprintf(1,'  ---  initialise model run \n\n');


%% setup numerical grid
% set dimensions of staggered/ghosted 2D grid
NUM.nxC         =  NUM.N+1;                             % number of corner nodes in x-direction
NUM.nzC         =  NUM.N+1;                             % number of corner nodes in z-direction
NUM.nxP         =  NUM.N+2;                            % number of centre nodes in x-direction
NUM.nzP         =  NUM.N+2;                            % number of centre nodes in z-direction
NUM.nxW         =  NUM.N+2;                            % number of z-face nodes in x-direction
NUM.nzW         =  NUM.N+1;                            % number of z-face nodes in z-direction
NUM.nxU         =  NUM.N+1;                            % number of x-face nodes in x-direction
NUM.nzU         =  NUM.N+2;                            % number of x-face nodes in z-direction
NUM.NC          =  NUM.nxC*NUM.nzC;                     % total number of corner nodes
NUM.NP          =  NUM.nxP*NUM.nzP;                     % total number of corner nodes
NUM.NW          =  NUM.nxW*NUM.nzW;                     % total number of z-face nodes
NUM.NU          =  NUM.nxU*NUM.nzU;                     % total number of x-face nodes
NUM.NDOF        =  NUM.NP+NUM.NW+NUM.NU;                % total number of all degrees of freedum

% set coordinate vectors
NUM.xC          =  0:NUM.h:NUM.L;                      % Horizontal coordinates of corner nodes [m]
NUM.zC          =  0:NUM.h:NUM.D;                      % Vertical   coordinates of corner nodes [m]
NUM.xP          =  -NUM.h/2:NUM.h:NUM.L+NUM.h/2;     % Horizontal coordinates of centre nodes [m]
NUM.zP          =  -NUM.h/2:NUM.h:NUM.D+NUM.h/2;     % Vertical   coordinates of centre nodes [m]
NUM.xW          =  NUM.xP;                              % Horizontal coordinates of z-face nodes [m]
NUM.zW          =  NUM.zC;                              % Vertical   coordinates of z-face nodes [m]
NUM.xU          =  NUM.xC;                              % Horizontal coordinates of x-face nodes [m]
NUM.zU          =  NUM.zP;                              % Vertical   coordinates of x-face nodes [m]

% set 2D coordinate grids
[NUM.XC,NUM.ZC] = meshgrid(NUM.xC,NUM.zC);              % corner nodes grid
[NUM.XP,NUM.ZP] = meshgrid(NUM.xP,NUM.zP);              % centre nodes grid
[NUM.XW,NUM.ZW] = meshgrid(NUM.xW,NUM.zW);              % z-face nodes grid
[NUM.XU,NUM.ZU] = meshgrid(NUM.xU,NUM.zU);              % x-face nodes grid

%% setup mapping arrays
NUM.Map         =  reshape(1:NUM.NP,NUM.nzP,NUM.nxP);
NUM.MapW        =  reshape(1:NUM.NW,NUM.nzW,NUM.nxW);
NUM.MapU        =  reshape(1:NUM.NU,NUM.nzU,NUM.nxU) + NUM.NW;
NUM.MapP        =  reshape(1:NUM.NP,NUM.nzP,NUM.nxP) + NUM.NW + NUM.NU;

% get smoothed initialisation field
rng(15);
rp = randn(NUM.nzP,NUM.nxP);
for i = 1:round(smth)
    rp(2:end-1,2:end-1) = rp(2:end-1,2:end-1) + diff(rp(:,2:end-1),2,1)./8 ...
                                              + diff(rp(2:end-1,:),2,2)./8;
    rp([1 end],:)       = rp([2 end-1],:);
    rp(:,[1 end])       = rp(:,[2 end-1]);
end
rp              = rp./max(abs(rp(:)));
rp              = rp - mean(mean(rp(2:end-1,2:end-1)));

%% setup material property arrays
% gravity
PHY.gxP = zeros(NUM.nzP,NUM.nxP) + PHY.gx;      PHY.gzP = zeros(NUM.nzP,NUM.nxP) + PHY.gz;
PHY.gx  = zeros(NUM.nzU,NUM.nxU) + PHY.gx;      PHY.gz  = zeros(NUM.nzW,NUM.nxW) + PHY.gz;
% materials
MAT.Eta         = zeros(NUM.nzP,NUM.nxP) + PHY.Eta0;   	% viscosity
MAT.Eta0        = MAT.Eta;                              % reference viscosity
MAT.aT          = zeros(NUM.nzP,NUM.nxP) + PHY.aT0;    	% thermal expansivity
MAT.kT          = zeros(NUM.nzP,NUM.nxP) + PHY.kT0;     % thermal conductivity
MAT.Cp          = zeros(NUM.nzP,NUM.nxP) + PHY.Cp0;   	% heat capacity
% MAT.Hr	 = zeros(NUM.nzP,NUM.nxP) + PHY.Hr0;	% radiogenic heating
MAT.Hr          = PHY.Hr0 + PHY.Hr0*0.01 .*rp;


%% set initial solution fields
SOL.phi         = PHY.x0 + (PHY.x1-PHY.x0) .* (1+erf(40*(NUM.ZP/NUM.D-zlay)))/2 + dc.*rp; %   meanT0 = mean(mean(T(2:end-1,2:end-1)));
SOL.phi(1,:)    = PHY.x0; SOL.phi(end,:) = PHY.x1;
MAT.Rhot       = SOL.phi.*PHY.RhoFe0 + (1-SOL.phi).*PHY.RhoSi0;

% temperature
pert = -NUM.h/2.*cos(NUM.XP*2*pi/NUM.D);
switch SOL.Ttype
    case 'constant'     % constant temperature
        SOL.T  = zeros(NUM.nzP,NUM.nxP) + SOL.T0;
    case 'linear'       % linear temperature gradient with depth
        SOL.T  = SOL.T0 + abs(NUM.ZP+pert)./NUM.D.*(SOL.T1-SOL.T0);
    case 'gaussian'     % constant temperature with gaussian plume
        SOL.T = zeros(NUM.nzP,NUM.nxP) + SOL.T0;
        SOL.T = SOL.T + SOL.dT.*exp(-(NUM.XP-SOL.xT).^2./SOL.rT.^2 - (NUM.ZP-SOL.zT).^2./SOL.rT.^2 );
    case 'hot bottom'
        SOL.T = zeros(NUM.nzP,NUM.nxP) + SOL.T0;
        SOL.T(end-10:end,:) = SOL.T0+100;
%         SOL.T = SOL.T + 1./(1+exp(-(NUM.ZP-SOL.zT+pert)./(NUM.D/50))) .* (SOL.T1-SOL.T0);
end
% convert from potential to natural temperature
SOL.T            = SOL.T.*exp(PHY.aT0*(PHY.gzP.*(NUM.ZP+pert) + PHY.gxP.*NUM.XP)./PHY.Cp0);
SOL.T([1 end],:) = SOL.T([2 end-1],:);
SOL.T(:,[1 end]) = SOL.T(:,[2 end-1]);

%% setup velocity-pressure solution arrays
SOL.W           = zeros(NUM.nzW,NUM.nxW);               % z-velocity on z-face nodes
SOL.U           = zeros(NUM.nzU,NUM.nxU);               % x-velocity on x-face nodes
SOL.P           = zeros(NUM.nzP,NUM.nxP);               % pressure on centre nodes
SOL.Stokes      = zeros(NUM.nzP,NUM.nxP);               % Stokes settling velocity
SOL.Pt          = MAT.Rhot.*(PHY.gzP.*NUM.ZP + PHY.gxP.*NUM.XP); % total pressure

% projected velocities on centre nodes
SOL.UP          = zeros(NUM.nzP,NUM.nxP);
SOL.WP          = zeros(NUM.nzP,NUM.nxP);
SOL.UP(:,2:end-1) = SOL.U(:,1:end-1)+SOL.U(:,2:end)./2;
SOL.WP(2:end-1,:) = SOL.W(1:end-1,:)+SOL.W(2:end,:)./2;
S               = [SOL.W(:);SOL.U(:);SOL.P(:)];             % Solution matrix

%% setup deformation property arrays
DEF.ups         = zeros(NUM.nzP,NUM.nxP);               % velocity divergence on centre nodes
DEF.exx         = zeros(NUM.nzP,NUM.nxP);               % x-normal strain rate on centre nodes
DEF.ezz         = zeros(NUM.nzP,NUM.nxP);               % z-normal strain rate on centre nodes
DEF.exz         = zeros(NUM.nzC,NUM.nxC);               % xz-shear strain rate on corner nodes
DEF.eII         = zeros(NUM.nzP,NUM.nxP);               % strain rate magnitude on centre nodes
DEF.txx         = zeros(NUM.nzP,NUM.nxP);               % x-normal stress on centre nodes
DEF.tzz         = zeros(NUM.nzP,NUM.nxP);               % z-normal stress on centre nodes
DEF.txz         = zeros(NUM.nzC,NUM.nxC);               % xz-shear stress on corner nodes
DEF.tII         = zeros(NUM.nzP,NUM.nxP);               % stress magnitude on centre nodes


%% setup heating rates
SOL.dTdt        = zeros(NUM.N  ,NUM.N );              % temperature rate of change
SOL.Hs          = zeros(NUM.nzP,NUM.nxP);          	% shear heating rate
SOL.Ha          = zeros(NUM.nzP,NUM.nxP);          	% adiabatic heating rate
Div_PHI         = zeros(NUM.N  ,NUM.N );              % level set free surface change


%% update nonlinear material properties
RUN.frame = 0;     % initialise output frame count
NUM.time  = 0;      % initialise time count
NUM.step  = 0;      % initialise time step count
up2date;


%% output initial condition
% output;

