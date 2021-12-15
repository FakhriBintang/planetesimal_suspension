% print update header
fprintf(1,'  ---  update materials & deformation \n');
tic;  % start clock on update

% update viscosity
MAT.EtaS                = (1-SOL.phi./0.3).^(-2);       % suspension viscosity
MAT.Mu                  = -4.13 +3703./(SOL.T-761.7);   % T-dependent melt viscosity
MAT.Eta                 = MAT.Eta0.*MAT.EtaS;           % mixture viscosity
MAT.EtaC                = (MAT.Eta(1:end-1,1:end-1) ...
                        +  MAT.Eta(2:end  ,1:end-1) ...
                        +  MAT.Eta(1:end-1,2:end  ) ...
                        +  MAT.Eta(2:end  ,2:end  ))/4;                     % interpolate to corner nodes

%% update T-dependent density
% silicate
MAT.RhoSi               = PHY.RhoSi0.*(1 - MAT.aT.*(SOL.T-SOL.T0));        % density on centre nodes
MAT.RhoUSi              = (MAT.RhoSi(:,1:end-1)+MAT.RhoSi(:,2:end))./2;   % density on the vx nodes
MAT.RhoWSi              = (MAT.RhoSi(1:end-1,:)+MAT.RhoSi(2:end,:))./2;   % density on the vz nodes

% iron droplets
MAT.RhoFe               = PHY.RhoFe0.*(1 - MAT.aT.*(SOL.T-SOL.T0));        % density on centre nodes
MAT.RhoUFe              = (MAT.RhoFe(:,1:end-1)+MAT.RhoFe(:,2:end))./2;   % density on the vx nodes
MAT.RhoWFe              = (MAT.RhoFe(1:end-1,:)+MAT.RhoFe(2:end,:))./2;   % density on the vz nodes

% Bulk density
MAT.Rhot                = SOL.phi.*MAT.RhoFe + (1-SOL.phi).*MAT.RhoSi;
MAT.RhoUt               = (MAT.Rhot(:,1:end-1)+MAT.Rhot(:,2:end))./2;     % density on the vx nodes
MAT.RhoWt               = (MAT.Rhot(1:end-1,:)+MAT.Rhot(2:end,:))./2;     % density on the vz nodes

% density contrast (between Fe and bulk)
MAT.dRho                = MAT.RhoFe - MAT.Rhot;

%% calculate stokes settling velocity  
if SOL.BCSeg==2; sds = -1;      % no slip
else;            sds = +1; end  % free slip

SOL.seg                 = 2.*MAT.dRho.*(PHY.d.^2).*PHY.gzP./9./MAT.Eta;
% set boundary conditions
SOL.seg(:,[1 end])      = sds*SOL.seg(:,[2 end-1]);
SOL.seg([1 end],:)      = 0;
SOL.segU = (SOL.seg(:,1:end-1)+SOL.seg(:,2:end))./2;
SOL.segW = (SOL.seg(1:end-1,:)+SOL.seg(2:end,:))./2;

if ~mod(NUM.step,round(2*RUN.nup/NUM.CFL))                          % only update deformation when fluid mechanics solved
    %% update physical time step
    VR                  = abs([SOL.U(:);SOL.W(:)]);                         % reference velocity
    dtadvR              = NUM.h/2/max(VR);
    VF                  = abs([SOL.U(:)+SOL.segU(:);SOL.W(:)+SOL.segW(:)]);              % Fe velocity
    dtadvF              = NUM.h/2/max(VF);                                  % 
    dtadvn              = min(dtadvR,dtadvF);                               % Stable timestep for advection
    
    kappa               = MAT.kT./MAT.Rhot./MAT.Cp;
    kappa               = kappa(:);
    dtdiff              = (NUM.h/2)^2 / max(kappa);            % stable time step for T diffusion
    
    NUM.dt              = NUM.CFL * min(dtdiff,dtadvn);                     % fraction of minimum stable time step
    if dtdiff<dtadvn 
    disp('diffusion regime')
    else
    disp('advection regime')
    end
    
    %% update strain-rate components
    % get volumetric strain-rate (velocity divergence)
    DEF.ups(2:end-1,2:end-1) = diff(SOL.U(2:end-1,:),1,2)./NUM.h ...
                             + diff(SOL.W(:,2:end-1),1,1)./NUM.h;           % velocity divergence
    DEF.ups([1 end],:)       = DEF.ups([2 end-1],:);                        % apply boundary conditions
    DEF.ups(:,[1 end])       = DEF.ups(:,[2 end-1]);
    
    % get deviatoric strain rates
    DEF.exx(:,2:end-1)      = diff(SOL.U,1,2)./NUM.h ...
                            - DEF.ups(:,2:end-1)./3;                        % x-normal strain rate
    DEF.exx([1 end],:)      = DEF.exx([2 end-1],:);                         % apply boundary conditions
    DEF.exx(:,[1 end])      = DEF.exx(:,[2 end-1]);
    
    DEF.ezz(2:end-1,:)      = diff(SOL.W,1,1)./NUM.h ...
                            - DEF.ups(2:end-1,:)./3;                        % z-normal strain rate
    DEF.ezz([1 end],:)      =  DEF.ezz([2 end-1],:);                        % apply boundary conditions
    DEF.ezz(:,[1 end])      =  DEF.ezz(:,[2 end-1]);
    
    DEF.exz                 = (diff(SOL.U,1,1)./NUM.h ...
                            +  diff(SOL.W,1,2)./NUM.h)/2;                   % shear strain rate
    
    % update strain-rate magnitude
    DEF.eII(2:end-1,2:end-1)= (  (DEF.exx(2:end-1,2:end-1).^2 + DEF.ezz(2:end-1,2:end-1).^2 ...
                            + 2.*(DEF.exz(1:end-1,1:end-1).^2 + DEF.exz(2:end,1:end-1).^2 ...
                            +     DEF.exz(1:end-1,2:end  ).^2 + DEF.exz(2:end,2:end  ).^2)/4 )/2).^0.5 + 1e-16;
    DEF.eII([1 end],:)      =     DEF.eII([2 end-1],:);                   	% apply boundaries
    DEF.eII(:,[1 end])      =     DEF.eII(:,[2 end-1]);
    
    
    %% update stress components
    DEF.txx = MAT.Eta  .* DEF.exx;                                          % x-normal stress
    DEF.tzz = MAT.Eta  .* DEF.ezz;                                          % z-normal stress
    DEF.txz = MAT.EtaC .* DEF.exz;                                          % xz-shear stress
    
    % update strain-rate magnitude
    DEF.tII(2:end-1,2:end-1)= (  (DEF.txx(2:end-1,2:end-1).^2 + DEF.tzz(2:end-1,2:end-1).^2 ...
                            + 2.*(DEF.txz(1:end-1,1:end-1).^2 + DEF.txz(2:end,1:end-1).^2 ...
                            +     DEF.txz(1:end-1,2:end  ).^2 + DEF.txz(2:end,2:end  ).^2)/4 )/2).^0.5 + 1e-16;
    DEF.tII([1 end],:)      =     DEF.tII([2 end-1],:);                     % apply boundaries
    DEF.tII(:,[1 end])      =     DEF.tII(:,[2 end-1]);
end
%% update heat source fields
    % update shear heating
    SOL.Hs = 2.*DEF.eII.*DEF.tII;
    
    % update adiabatic heating
    SOL.Ha = ((SOL.WP.*PHY.gzP + SOL.UP.*PHY.gxP) .* MAT.Rhot + SOL.seg.*PHY.gzP).*MAT.aT.*SOL.T;