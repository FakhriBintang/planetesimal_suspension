clear all; close all
RunID = 'bnchm_VP_mms';

[~,systemname]  = system('hostname');
systemname(end) = [];
switch systemname
    case 'Horatio'
        outpath = ['/media/43TB_RAID_Array/fbintang/test_out/out/',RunID];
        if ~exist(outpath, 'dir'); mkdir(outpath); end
    otherwise
        outpath = ['../out/',RunID];
        if ~exist(outpath, 'dir'); mkdir(outpath); end
end

L = 10; D = 10;
g0 = 0.1;
NN = [50, 100, 200];
for N = NN
    nop = 1;
    bnchm = 1;

    h = D/N; % grid spacing

    %inner indeces
    inx = 2:N-1;
    inz = 2:N-1;

    % create manufactured solution
    clear x z
    TINY = 1e-16;
    syms U_mms(x,z) W_mms(x,z) P_mms(x,z) eta_mms(x,z) rho_mms(x,z) src_mms(x,z)

    % compose manufactured solution variables
    W_mms(x,z) = 5.00e-5.*(cos(2*pi/(L/2)*(x)).*sin(2*pi/(L/2)*(z)));
    U_mms(x,z) = 1.00e-5.*(sin(2*pi/(L/2)*(x)).*cos(2*pi/(L/2)*(z)));
    P_mms(x,z) = 2.00e+3.*(sin(2*pi/(L/2)*(x)).*sin(2*pi/(L/2)*(z)));

    % compose manufactured material coefficients and volume source
    eta_mms(x,z)  = 1e+3-8e+2.*(cos(4*(x)*pi/L).*sin(4*(z)*pi/L));
    rho_mms(x,z)  = 3e+3-1e+1.*(cos(4*(x)*pi/L).*sin(4*(z)*pi/L)); rhoref = 3e+3;
    src_mms(x,z)  =     -1e-5.*(sin(4*(x)*pi/L).*sin(4*(z)*pi/L));

    % update strain rates
    DivV_mms(x,z)= (diff(W_mms,z) + diff(U_mms,x));
    exx_mms(x,z) = diff(U_mms,x) - DivV_mms./2;         % x-normal strain rate
    ezz_mms(x,z) = diff(W_mms,z) - DivV_mms./2;         % z-normal strain rate
    exz_mms(x,z) = 1/2.*(diff(U_mms,z)+diff(W_mms,x));  % xz-shear strain rate
    fprintf(1,' . ');

    % update stresses
    txx_mms(x,z) = eta_mms .* exx_mms;                  % x-normal stress
    tzz_mms(x,z) = eta_mms .* ezz_mms;                  % z-normal stress
    txz_mms(x,z) = eta_mms .* exz_mms;                  % xz-shear stress
    fprintf(1,' . ');

    % manufactured solution residuals
    res_W_mms = (diff(tzz_mms,z) + diff(txz_mms,x)) - diff(P_mms,z) + (rho_mms(x,z)-rhoref)*g0;
    res_U_mms = (diff(txx_mms,x) + diff(txz_mms,z)) - diff(P_mms,x);
    res_P_mms =-(diff(  W_mms,z) + diff(  U_mms,x))                 + src_mms(x,z);
    fprintf(1,' . ');

    % plot manufactured solution
    figure(15);
%     colormap(ocean);
    subplot(2,3,1); fcontour( -W_mms  ,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufcat. $W$ [m/hr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    subplot(2,3,2); fcontour(  U_mms  ,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufact. $U$ [m/hr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    subplot(2,3,3); fcontour(  P_mms/1e3 ,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufact. $P$ [kPa]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    fprintf(1,' . ');
    subplot(2,3,4); fcontour(      rho_mms ,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufact. $\rho$ [kg/m$^3$]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    subplot(2,3,5); fcontour(log10(eta_mms),[0,L],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufact. $\eta$ [log$_{10}$ Pas]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    subplot(2,3,6); fcontour(      src_mms ,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufact. $\dot{V} [1/s]$','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    drawnow;
    fprintf(1,' . \n');

    % evaluate mms source terms on appropriate coordinate grids
    fprintf(1,'\n  ***  evaluate manufactured solution\n\n');
    x_mms  = -h/2:h:L+h/2;
    z_mms  = -h/2:h:L+h/2;
    xu_mms = (x_mms(1:end-1)+x_mms(2:end))./2;
    zw_mms = (z_mms(1:end-1)+z_mms(2:end))./2;

    fprintf(1,'       Patience, my young Padawan!\n');
    fprintf(1,'       . ');

    [x,z] = meshgrid(x_mms,zw_mms);
    src_W_mms = double(subs(res_W_mms)); fprintf(1,' . ');
    [x,z] = meshgrid(xu_mms,z_mms);
    src_U_mms = double(subs(res_U_mms)); fprintf(1,' . ');
    [x,z] = meshgrid(x_mms,z_mms);
    src_P_mms = double(subs(res_P_mms)); fprintf(1,' . ');

    % plot manufactured residuals and evaluated source terms
    figure(16);
%     colormap(ocean);
    subplot(2,3,1); fcontour(-res_W_mms,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar; box on; title('manufactured $W$-res','Interpreter','latex');
    subplot(2,3,2); fcontour(-res_U_mms,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar; box on; title('manufactured $U$-res','Interpreter','latex');
    subplot(2,3,3); fcontour(-res_P_mms,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar; box on; title('manufactured $P$-res','Interpreter','latex');
    subplot(2,3,4); imagesc(x_mms,zw_mms,-src_W_mms); axis ij equal tight; colorbar; box on; title('evaluated $W$-res','Interpreter','latex');
    subplot(2,3,5); imagesc(xu_mms,z_mms,-src_U_mms); axis ij equal tight; colorbar; box on; title('evaluated $U$-res','Interpreter','latex');
    subplot(2,3,6); imagesc(x_mms ,z_mms,-src_P_mms); axis ij equal tight; colorbar; box on; title('evaluated $P$-res','Interpreter','latex');
    drawnow;

    % evaluate analytical solution on appropriate coordinate grids
    [x,z]  = meshgrid(x_mms,zw_mms);
    W_mms  = double(subs(W_mms)); fprintf(1,' . ');
    rhoBF  = double(subs(rho_mms)); fprintf(1,' . ');
    rhoBF  = rhoBF(2:end-1,2:end-1);
    [x,z]  = meshgrid(xu_mms,z_mms);
    U_mms  = double(subs(U_mms)); fprintf(1,' . ');
    [x,z]  = meshgrid(x_mms,z_mms);
    P_mms  = double(subs(P_mms)); fprintf(1,' . ');
    Eta    = double(subs(eta_mms)); fprintf(1,' . ');
%     Eta    = Eta(2:end-1,2:end-1);
%     gz = zeros(Nz+1,Nx+2) +g0;
    VolSrc = double(subs(src_mms)); fprintf(1,' . ');
    VolSrc = VolSrc(2:end-1,2:end-1);
    [x,z]  = meshgrid(xu_mms,zw_mms);
    EtaC  = double(subs(eta_mms)); fprintf(1,' . ');

    WBG    = 0.*W_mms;
    UBG    = 0.*U_mms;
    SOL    = [W_mms(:);U_mms(:);P_mms(:)];

    % get mapping arrays
    Nx = length(x_mms)-2;
    Nz = length(z_mms)-2;
    nzP = Nz+2; nxP = Nx+2;
    nzW = Nz+1; nxW = Nx+2;
    nzU = Nz+2; nxU = Nx+1;
    NP = nzP * nxP;
    NW = nzW * nxW;
    NU = nzU * nxU;
    MapP = reshape(1:NP,Nz+2,Nx+2);
    MapW = reshape(1:NW,Nz+1,Nx+2);
    MapU = reshape(1:NU,Nz+2,Nx+1) + NW;
    
    gz = zeros(Nz+1,Nx+2) +g0;
    % set boundary conditions to free slip
    BCsides = -1;
    BCtop = -1;
    BCbot = -1;

    fprintf(1,' . \n');

    % call fluid mechanics solver
    FMtime = 0;
    run('../src/solve_fluidmech');

    figure(17); clf;
%     colormap(ocean);
    subplot(2,3,1); imagesc(x_mms,zw_mms,-W(:,inx)); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. sol. $W$ [m/hr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    subplot(2,3,2); imagesc(xu_mms,z_mms, U(inz,:)); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. sol. $U$ [m/hr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    subplot(2,3,3); imagesc(x_mms ,z_mms, P(inz,inx)/1e3); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. sol. $P$ [kPa]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    subplot(2,3,4); imagesc(x_mms,zw_mms,-(W(:,inx)-W_mms(:,inx))); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. err. $W$ [m/hr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    subplot(2,3,5); imagesc(xu_mms,z_mms, (U(inz,:)-U_mms(inz,:))); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. err. $U$ [m/hr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    subplot(2,3,6); imagesc(x_mms ,z_mms, (P(inz,inx)-P_mms(inz,inx))/1e3); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. err. $P$ [kPa]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    drawnow;

    % get solution error
    EW = norm(W(:  ,inx)-W_mms(:  ,inx),'fro')./norm(W_mms(:  ,inx),'fro');
    EU = norm(U(inz,:  )-U_mms(inz,:  ),'fro')./norm(U_mms(inz,:  ),'fro');
    EP = norm(P(inz,inx)-P_mms(inz,inx),'fro')./norm(P_mms(inz,inx),'fro');
    
    % plot error convergence
    fh18 = figure(18);
    p1 = loglog(h,EW,'rs','MarkerSize',8,'LineWidth',2); axis xy tight; hold on; box on;
    p2 = loglog(h,EU,'go','MarkerSize',8,'LineWidth',2); axis xy tight; hold on; box on;
    p3 = loglog(h,EP,'bv','MarkerSize',8,'LineWidth',2); axis xy tight; hold on; box on;
    set(gca,'TicklabelInterpreter','latex','FontSize',12)
    xlabel('grid step [m]','Interpreter','latex')
    ylabel('rel. numerical error [1]','Interpreter','latex')
    set(gca,'TicklabelInterpreter','latex')
    title('Numerical convergence in space','Interpreter','latex','FontSize',20)

    if N == NN(1)
        p4 = loglog(D./NN,mean([EW,EU,EP]).*(NN(1)./NN).^2,'k-','LineWidth',2);  % plot linear trend for comparison
    end
    if N == NN(end)
        legend([p1,p2,p3,p4],{'error W','error U','error P','quadratic'},'Interpreter','latex','box','on','location','southeast')
    end

%     % plot error convergence
%     fh19 = figure(19);
%     DOFS = (NN+2).*(NN+2) + 2.*(NN+1).*(NN+2);
%     dofs = (N+2).*(N+2) + 2.*(N+1).*(N+2);
%     p5 = loglog(dofs,FMtime,'r+','MarkerSize',8,'LineWidth',2); axis xy tight; hold on; box on;
%     set(gca,'TicklabelInterpreter','latex','FontSize',12)
%     xlabel('\# dofs [1]','Interpreter','latex','FontSize',16)
%     ylabel('time to solution [s]','Interpreter','latex','FontSize',16)
%     title('Scaling of direct solver','Interpreter','latex','FontSize',20)
% 
%     if N == NN(1)
%         p6 = loglog(DOFS,0.95*FMtime*(DOFS./DOFS(1)).^1,'k-','LineWidth',2);  % plot linear trend for comparison
%     end
%     if N == NN(end)
%         legend([p5,p6],{'time to solution','linear'},'Interpreter','latex','box','on','location','southeast')
%     end   

end
name = [outpath,'/',RunID,'_bnchm'];
print(fh18,name,'-dpng','-r300','-vector');