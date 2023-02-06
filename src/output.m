% planetesimal: generate output

% print output header
fprintf(1,'\n*****  prepare output frame %d  *****\n',step/nop);


%% plot output figures
if plot_op
    
    % prepare for plotting
    TX = {'Interpreter','Latex'}; FS = {'FontSize',12};
    TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',10};
    UN = {'Units','Centimeters'};
    
    if Nx <= 10 && Nz <= 10  % create 0D plots
        
        fh1 = figure(1); clf;
        subplot(4,1,1)
        plot(HST.time/yr,HST.T(:,2),'LineWidth',2); axis xy tight; box on;
        title('$T [^\circ$C]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(4,1,2)
        plot(HST.time/yr,HST.xFe(:,2).*100,'LineWidth',2); axis xy tight; box on;
        title('$x_{Fe}$ [wt\% Fe-FeS]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(4,1,3)
        plot(HST.time/yr,HST.cFe(:,2).*100,'LineWidth',2); axis xy tight; box on;
        title('$\bar{c}_{Fe}$ [wt\% S]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(4,1,4)
        plot(HST.time/yr,HST.cSi(:,2).*100,'LineWidth',2); axis xy tight; box on;
        title('$\bar{c}_{Si}$ [wt\% SiO$_2$]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        xlabel('Time [yr]',TX{:},FS{:});
        
        fh2 = figure(2); clf;
        subplot(4,1,1)
        plot(HST.time/yr,HST.flFe(:,2).*100,'LineWidth',2); axis xy tight; box on;
        title(['$f_{Fe}$ [wt\%]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(4,1,2)
        plot(HST.time/yr,HST.flSi(:,2).*100,'LineWidth',2); axis xy tight; box on;
        title(['$f_{Si}$ [wt\%]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(4,1,3)
        plot(HST.time/yr,HST.rho(:,2),'LineWidth',2); axis xy tight; box on;
        title('$\bar{\rho}$ [kg/m$^3$]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(4,1,4)
        plot(HST.time/yr,log10(HST.Eta(:,2)),'LineWidth',2); axis xy tight; box on;
        title('$\bar{\Eta}$ [log$_{10}$ Pas]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        xlabel('Time [yr]',TX{:},FS{:});
        
    elseif Nx <= 10 % 1D plots
        
        fh1 = figure(1); clf;
        subplot(1,4,1)
        plot(mean(T(2:end-1,2:end-1),2),zP(2:end-1),'LineWidth',2); axis ij tight; box on;
        title('$T [^\circ$C]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        ylabel('Depth [m]',TX{:},FS{:});
        subplot(1,4,2)
        plot(mean(xFe(2:end-1,2:end-1),2).*100,zP(2:end-1),'LineWidth',2); axis ij tight; box on;
        title('$x_{Fe}$ [wt\% Fe-FeS]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,4,3)
        plot(mean(clFe(2:end-1,2:end-1),2).*100,zP(2:end-1),'-r','LineWidth',2); axis ij tight; box on; hold on
        plot(mean(csFe(2:end-1,2:end-1),2).*100,zP(2:end-1),'-b','LineWidth',2);
        plot(mean(cFe(2:end-1,2:end-1),2).*100,zP(2:end-1),'k-','LineWidth',2);
        title('$\bar{c}_{Fe}$ [wt\% S]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,4,4)
        plot(mean(clSi(2:end-1,2:end-1),2).*100,zP(2:end-1),'-r','LineWidth',2); axis ij tight; box on; hold on
        plot(mean(csSi(2:end-1,2:end-1),2).*100,zP(2:end-1),'-b','LineWidth',2);
        plot(mean(cSi(2:end-1,2:end-1),2).*100,zP(2:end-1),'k-','LineWidth',2);
        title('$\bar{c}_{Si}$ [wt\% SiO$_2$]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        
        fh2 = figure(2); clf;
        subplot(1,4,1)
        plot(mean(flFe(2:end-1,2:end-1),2).*100,zP(2:end-1),'LineWidth',2); axis ij tight; box on; hold on
        plot(mean(flSi(2:end-1,2:end-1),2).*100,zP(2:end-1),'LineWidth',2);
        title(['$f_{j}^\ell$ [wt\%]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
        legend('Fe', 'Si')
        ylabel('Depth [m]',TX{:},FS{:});
        subplot(1,4,2)
        plot(mean(phisFe(2:end-1,2:end-1),2).*100,zP(2:end-1),'LineWidth',2); axis ij tight; box on; hold on
        plot(mean(philFe(2:end-1,2:end-1),2).*100,zP(2:end-1),'LineWidth',2);
        plot(mean(phisSi(2:end-1,2:end-1),2).*100,zP(2:end-1),'LineWidth',2);
        plot(mean(philSi(2:end-1,2:end-1),2).*100,zP(2:end-1),'LineWidth',2);
        title(['$\phi_{j}^i$ [vol\%]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
        legend('\phi_{Fe}^s','\phi_{Fe}^{l}','\phi_{Si}^s','\phi_{Si}^{l}')
        subplot(1,4,3)
        plot(mean(GFe(2:end-1,2:end-1)./rho(2:end-1,2:end-1).*100.*yr,2),zP(2:end-1),'LineWidth',2); axis ij tight; box on;
        title('$\Gamma_{Fe}$ [wt\%/yr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,4,4)
        plot(mean(GSi(2:end-1,2:end-1)./rho(2:end-1,2:end-1).*100.*yr,2),zP(2:end-1),'LineWidth',2); axis ij tight; box on;
        title('$\Gamma_{Si}$ [wt\%/yr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        
        fh3 = figure(3); clf;
        subplot(1,4,1)
        plot(mean(-(philFe(1:end-1,2:end-1)+philFe(2:end,2:end-1))/2.*seglFe(:,2:end-1),2)*yr,zW,'LineWidth',2); hold on
        plot(mean(-(phisFe(1:end-1,2:end-1)+phisFe(2:end,2:end-1))/2.*segsFe(:,2:end-1),2)*yr,zW,'LineWidth',2);
        plot(mean(-(philSi(1:end-1,2:end-1)+philSi(2:end,2:end-1))/2.*seglSi(:,2:end-1),2)*yr,zW,'LineWidth',2);
        plot(mean(-(phisSi(1:end-1,2:end-1)+phisSi(2:end,2:end-1))/2.*segsSi(:,2:end-1),2)*yr,zW,'LineWidth',2);
        plot(mean(-W(:,2:end-1),2)*yr,zW,'k','LineWidth',2); axis ij tight; box on;
        title('$W$ [m/yr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        ylabel('Depth [m]',TX{:},FS{:});
        subplot(1,4,2)
        plot(mean(P(2:end-1,2:end-1),2),zP(2:end-1),'LineWidth',2); axis ij tight; box on;
        title('$P$ [Pa]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,4,3)
        plot(mean(rho(2:end-1,2:end-1),2),zP(2:end-1),'LineWidth',2); axis ij tight; box on;
        title('$\bar{\rho}$ [kg/m$^3$]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,4,4)
        plot(mean(log10(Eta(2:end-1,2:end-1)),2),zP(2:end-1),'LineWidth',2); axis ij tight; box on;
        title('$\bar{\eta}$ [log$_{10}$ Pas]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        
        fh4 = figure(4); clf;
        subplot(1,4,1)
        plot(mean(-(philFe(1:end-1,2:end-1)+philFe(2:end,2:end-1))/2.*seglFe(:,2:end-1),2)*yr,zW,'LineWidth',2); axis ij tight; box on;
        title(['$w_{\Delta,Fe}^\ell$ [m/yr]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
        ylabel('Depth [m]',TX{:},FS{:});
        subplot(1,4,2)
        plot(mean(-(phisFe(1:end-1,2:end-1)+phisFe(2:end,2:end-1))/2.*segsFe(:,2:end-1),2)*yr,zW,'LineWidth',2); axis ij tight; box on;
        title(['$w_{\Delta,Fe}^s$ [m/yr]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,4,3)
        plot(mean(-(philSi(1:end-1,2:end-1)+philSi(2:end,2:end-1))/2.*seglSi(:,2:end-1),2)*yr,zW,'LineWidth',2); axis ij tight; box on;
        title(['$w_{\Delta,Si}^\ell$ [m/yr]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,4,4)
        plot(mean(-(phisSi(1:end-1,2:end-1)+phisSi(2:end,2:end-1))/2.*segsSi(:,2:end-1),2)*yr,zW,'LineWidth',2); axis ij tight; box on;
        title(['$w_{\Delta,Si}^s$ [m/yr]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
        
        % plot conserved quantities
        fh5 = figure(5); clf;
        subplot(1,4,1)
        plot(mean(XFe(2:end-1,2:end-1),2),zP(2:end-1),'LineWidth',2); axis ij tight; box on; hold on;
        plot(mean(XSi(2:end-1,2:end-1),2),zP(2:end-1),'LineWidth',2);
        plot(mean(rho(2:end-1,2:end-1),2),zP(2:end-1),'-k','LineWidth',2);
        plot(mean(XFe(2:end-1,2:end-1),2)+ mean(XSi(2:end-1,2:end-1),2),zP(2:end-1),'--','LineWidth',2);
        legend('XFe', 'XSi', 'rho', 'XFe +XSi')
        ylabel('Depth [m]',TX{:},FS{:});
        subplot(1,4,2)
        plot(mean(CFe(2:end-1,2:end-1),2),zP(2:end-1),'LineWidth',2); axis ij tight; box on;
        title('$C_{Fe}$ ',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,4,3)
        plot(mean(CSi(2:end-1,2:end-1),2),zP(2:end-1),'LineWidth',2); axis ij tight; box on;
        title('$C_{Si}$ ',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,4,4)
        plot(mean(FlFe(2:end-1,2:end-1),2),zP(2:end-1),'LineWidth',2); axis ij tight; box on; hold on;
        plot(mean(FsFe(2:end-1,2:end-1),2),zP(2:end-1),'LineWidth',2);
        plot(mean(FlSi(2:end-1,2:end-1),2),zP(2:end-1),'LineWidth',2);
        plot(mean(FsSi(2:end-1,2:end-1),2),zP(2:end-1),'LineWidth',2);
        legend('FlFe','FsFe' , 'FlSi', 'FsSi')
        ylabel('Depth [m]',TX{:},FS{:});
        
    else
        
        % prepare for plotting
        TX = {'Interpreter','Latex'}; FS = {'FontSize',12};
        TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',10};
        UN = {'Units','Centimeters'};
        xq = round(N/10/2)+1:round(N/10):nxP;  % x-indexes for quiver plots
        zq = round(N/10/2)+1:round(N/10):nxP;  % z-indexes for quiver plots
        
        
        %% plot velocity and temperature solutions
        fh1  = figure(1); clf
        figure(1)
        subplot(2,3,1)
        imagesc(xU(:),zU(2:end-1),U(2:end-1,:)*yr); hold on;
        quiver(xP(xq),zP(zq),UP(zq,xq),WP(zq,xq),'k')
        colormap(subplot(2,3,1),cm2)
        axis ij equal tight;
        colorbar
        title('v^*_x [m/yr]')
        
        subplot(2,3,2)
        imagesc(xW(2:end-1),zW(:),-W(:,2:end-1)*yr); hold on;
        quiver(xP(xq),zP(zq),UP(zq,xq),WP(zq,xq),'k')
        colormap(subplot(2,3,2),cm2)
        axis ij equal tight;
        colorbar
        title('v_z^*-velocity [m/yr]')
        
        subplot(2,3,3)
        imagesc(xP(2:end-1),zP(2:end-1),P(2:end-1,2:end-1)); hold on;
        colormap(subplot(2,3,3),cm2)
        axis ij equal tight;
        colorbar
        title('dynamic pressure [Pa]')
        
        subplot(2,3,4)
        imagesc(xP(2:end-1),zP(2:end-1),T(2:end-1,2:end-1));
        colormap(subplot(2,3,4),flipud(cm1))
        axis ij equal tight;
        colorbar
        title('Temperature [C]')
        
        subplot(2,3,5);
        imagesc(xP(2:end-1),zP(2:end-1),rho(2:end-1,2:end-1));
        colormap(subplot(2,3,5),cm1)
        axis ij equal tight;
        colorbar
        title('bulk density [kgm^-^3]')
        
        subplot(2,3,6);
        imagesc(xP(2:end-1),zP(2:end-1),log10(Eta(2:end-1,2:end-1)));
        colormap(subplot(2,3,6),cm1)
        axis ij equal tight;
        colorbar
        title('\Eta [log_1_0 Pas]')
        
        
        %% plot phase vol fractions
        fh2 = figure(2); clf;
        figure(2) % plot phase fractions
        subplot(2,2,1);
        imagesc(xP(2:end-1),zP(2:end-1),philSi(2:end-1,2:end-1));
        colorbar
        axis ij equal tight;
        title('\phi_{Si}^l')
        subplot(2,2,2);
        imagesc(xP(2:end-1),zP(2:end-1),phisSi(2:end-1,2:end-1));
        colorbar
        axis ij equal tight;
        title('\phi_{Si}^s')
        subplot(2,2,3);
        imagesc(xP(2:end-1),zP(2:end-1),philFe(2:end-1,2:end-1));
        colorbar
        axis ij equal tight;
        title('\phi_{Fe}^l')
        subplot(2,2,4);
        imagesc(xP(2:end-1),zP(2:end-1),phisFe(2:end-1,2:end-1));
        colorbar
        axis ij equal tight;
        title('\phi_{Fe}^s')
        
        
        %% plot phase
        fh3 = figure(3); clf
        figure(3) % plot phase fractions
        subplot(2,2,1);
        imagesc(xP(2:end-1),zP(2:end-1),flSi(2:end-1,2:end-1));
        colorbar
        axis ij equal tight
        title('f_{Si}^l [wt]')
        subplot(2,2,2);
        imagesc(xP(2:end-1),zP(2:end-1),fsSi(2:end-1,2:end-1));
        colorbar
        axis ij equal tight
        title('f_{Si}^s [wt]')
        subplot(2,2,3);
        imagesc(xP(2:end-1),zP(2:end-1),flFe(2:end-1,2:end-1));
        colorbar
        axis ij equal tight
        title('f_{Fe}^l [wt]')
        subplot(2,2,4);
        imagesc(xP(2:end-1),zP(2:end-1),fsFe(2:end-1,2:end-1));
        colorbar
        axis ij equal tight
        title('f_{Fe}^s [wt]')
        
        
        %% plot phase segregation
        fh4 = figure(4); clf
        figure(4) % plot phase segregation
        subplot(2,2,1);
        imagesc(xW(2:end-1),zW(:),-seglSi(:,2:end-1)*yr);
        colorbar
        axis ij equal tight
        title('\Delta v_{Si}^l [m/yr]')
        
        subplot(2,2,2);
        imagesc(xW(2:end-1),zW(:),-segsSi(:,2:end-1)*yr);
        colorbar
        axis ij equal tight
        title('\Delta v_{Si}^s [m/yr]')
        
        subplot(2,2,3);
        imagesc(xW(2:end-1),zW(:),-seglFe(:,2:end-1)*yr);
        colorbar
        axis ij equal tight
        title('\Delta v_{Fe}^l [m/yr]')
        
        subplot(2,2,4);
        imagesc(xW(2:end-1),zW(:),-segsFe(:,2:end-1)*yr);
        colorbar
        axis ij equal tight
        title('\Delta v_{Fe}^s [m/yr]')
        
        
        %% plot system fractions and chemical solutions
        fh5 = figure(5); clf
        figure(5)
        subplot(2,2,1)
        imagesc(xP(2:end-1),zP(2:end-1),xFe(2:end-1,2:end-1))
        colorbar
        axis ij equal tight
        title('x_{Fe} [wt]')
        
        subplot(2,2,2)
        imagesc(xP(2:end-1),zP(2:end-1),xSi(2:end-1,2:end-1))
        colorbar
        axis ij equal tight
        title('x_{Si} [wt]')
        
        subplot(2,2,3)
        imagesc(xP(2:end-1),zP(2:end-1),cFe(2:end-1,2:end-1))
        colorbar
        axis ij equal tight
        title('c_{Fe} [wt]')
        
        subplot(2,2,4)
        imagesc(xP(2:end-1),zP(2:end-1),cSi(2:end-1,2:end-1))
        colorbar
        axis ij equal tight
        title('c_{Si} [wt]')
        
        
        %% plot phase densities
        fh6 = figure(6); clf
        figure(6)
        subplot(2,2,1)
        imagesc(xP(2:end-1),zP(2:end-1),rholSi(2:end-1,2:end-1))
        colorbar
        axis ij equal tight
        title('\rho_{Si}^l')
        
        subplot(2,2,2)
        imagesc(xP(2:end-1),zP(2:end-1),rhosSi(2:end-1,2:end-1))
        colorbar
        axis ij equal tight
        title('\rho_{Si}^s')
        
        subplot(2,2,3)
        imagesc(xP(2:end-1),zP(2:end-1),rholFe(2:end-1,2:end-1))
        colorbar
        axis ij equal tight
        title('\rho_{Fe}^l')
        
        subplot(2,2,4)
        imagesc(xP(2:end-1),zP(2:end-1),rhosFe(2:end-1,2:end-1))
        colorbar
        axis ij equal tight
        title('\rho_{Fe}^s')
        
        
        %% plot chemical phase compositions
        fh7 = figure(7); clf
        figure(7)
        subplot(2,2,1)
        imagesc(xP(2:end-1),zP(2:end-1),clSi(2:end-1,2:end-1))
        colorbar
        axis ij equal tight
        title('c_{Si}^l')
        
        subplot(2,2,2)
        imagesc(xP(2:end-1),zP(2:end-1),csSi(2:end-1,2:end-1))
        colorbar
        axis ij equal tight
        title('c_{Si}^s')
        
        subplot(2,2,3)
        imagesc(xP(2:end-1),zP(2:end-1),clFe(2:end-1,2:end-1))
        colorbar
        axis ij equal tight
        title('c_{Fe}^l')
        
        subplot(2,2,4)
        imagesc(xP(2:end-1),zP(2:end-1),csFe(2:end-1,2:end-1))
        colorbar
        axis ij equal tight
        title('c_{Fe}^s')
        
        
%         %% plot Partial densities and
%         fh8 = figure(8); clf;
%         figure(8)
%         subplot(2,2,1)
%         imagesc(xP(2:end-1),zP(2:end-1),CFe(2:end-1,2:end-1))
%         colorbar
%         axis ij equal tight
%         title('C_{Fe} ')
%         
%         subplot(2,2,2)
%         imagesc(xP(2:end-1),zP(2:end-1),CSi(2:end-1,2:end-1))
%         colorbar
%         axis ij equal tight
%         title('C_{Si} ')
%         
%         subplot(2,2,3)
%         imagesc(xP(2:end-1),zP(2:end-1),XFe(2:end-1,2:end-1))
%         colorbar
%         axis ij equal tight
%         title('X_{Fe}')
%         
%         subplot(2,2,4)
%         imagesc(xP(2:end-1),zP(2:end-1),XSi(2:end-1,2:end-1))
%         colorbar
%         axis ij equal tight
%         title('X_{Si}')
        
    end
    
    %% plot phase diagrams
    fh9 = figure(9); clf;
    TT = linspace(TSi1,TSi2,1e3);
    cc = [linspace(cphsSi2,(perCsSi+perClSi)/2,round((perTSi-TSi1)./(TSi2-TSi1)*1e3)),linspace((perCsSi+perClSi)/2,cphsSi1,round((perTSi-TSi2)./(TSi1-TSi2)*1e3))];
    [~,CCxSi,CClSi] = equilibrium(TT,cc,0.*TT,TSi1,TSi2,cphsSi1,cphsSi2,...
                                  perTSi,perCsSi,perClSi,clap,PhDgSi);
    
    TT2 = linspace(TFe1,TFe2,1000);
    cc2 = linspace(cphsFe2,cphsFe1,length(TT2));
    [~,CCxFe,CClFe] = equilibrium(TT2,cc2,0.*TT2,TFe1,TFe2,cphsFe1,cphsFe2,...
                                  perTFe,perCsFe,perClFe,clap,PhDgFe);
    subplot(1,2,1)
    plot(CCxSi,TT,'k-','LineWidth',2); axis tight; axis square; hold on; box on;
    plot(CClSi,TT,'k-','LineWidth',2); axis tight; hold on; axis square; box on;
    plot(cSi(2:end-1,2:end-1),T(2:end-1,2:end-1), '.k', csSi(2:end-1,2:end-1),T(2:end-1,2:end-1), '.b', clSi(2:end-1,2:end-1),T(2:end-1,2:end-1), '.r','MarkerSize',25)
    %title('SiO$_2$ Phase Diagram','Interpreter','latex','FontSize',18)
    xlabel('Major component [wt\% SiO$_2$]','Interpreter','latex','FontSize',15)
    ylabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)

    subplot(1,2,2)
    plot(CCxFe,TT2,'k-','LineWidth',2); axis tight; axis square; hold on; box on;
    plot(CClFe,TT2,'k-','LineWidth',2); axis tight; axis square; hold on; box on;
    plot(cFe(2:end-1,2:end-1),T(2:end-1,2:end-1), '.k', csFe(2:end-1,2:end-1),T(2:end-1,2:end-1), '.b', clFe(2:end-1,2:end-1),T(2:end-1,2:end-1), '.r','MarkerSize',25)
    %title('Fe-FeS Phase Diagram','Interpreter','latex','FontSize',18)
    xlabel('Major component [wt\% S]','Interpreter','latex','FontSize',15)
    ylabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
    drawnow;
    
    
    %% plot conserved quantities
    fh10 = figure(10); clf;
    % mass/energy conservation
    subplot(5,1,1)
    plot(HST.time./yr,HST.EM  ,'k-','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
    ylabel('consv. $Mass$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(5,1,2)
    plot(HST.time./yr,HST.ES  ,'k-','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
    ylabel('consv. $S$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(5,1,3)
    plot(HST.time./yr,HST.EXFe,'k-','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
    ylabel('consv. $X_{Fe}$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(5,1,4)
    plot(HST.time./yr,HST.ECFe,'k-','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
    ylabel('consv. $C_{Fe}$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(5,1,5)
    plot(HST.time./yr,HST.ECSi,'k-','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
    ylabel('consv. $C_{Si}$',TX{:},FS{:}); set(gca,TL{:},TS{:})
    xlabel('Time [yr]',TX{:},FS{:});
    
end


%% save output
if save_op
    if Nx <= 10 && Nz <= 10  % save 0D plots
        
        name = [outpath '/',RunID,'_txc_',num2str(floor(step/nop))]; % figure 1
        print(fh1,name,'-dpng','-r300','-opengl');
        name = [outpath '/',RunID,'_fre_',num2str(floor(step/nop))]; % figure 2
        print(fh2,name,'-dpng','-r300','-opengl');
        
    elseif Nx <= 10   % save 1D plots
        
        name = [outpath '/',RunID,'_txc_',num2str(floor(step/nop))];   % figure 1
        print(fh1,name,'-dpng','-r300','-opengl');
        name = [outpath '/',RunID,'_fmr_',num2str(floor(step/nop))];   % figure 2
        print(fh2,name,'-dpng','-r300','-opengl');
        name = [outpath '/',RunID,'_vp_',num2str(floor(step/nop))];    % figure 3
        print(fh3,name,'-dpng','-r300','-opengl');
        name = [outpath '/',RunID,'_seg_',num2str(floor(step/nop))];   % figure 4
        print(fh4,name,'-dpng','-r300','-opengl');
        name = [outpath '/',RunID,'_consv_',num2str(floor(step/nop))]; % figure 5
        print(fh5,name,'-dpng','-r300','-opengl');
        
    else  % save 2D plots
        
        name = [outpath '/',RunID,'_vp_',num2str(floor(step/nop))]; % figure 1
        print(fh1,name,'-dpng','-r300','-opengl');
        name = [outpath '/',RunID,'_phsvol_',num2str(floor(step/nop))]; % figure 2
        print(fh2,name,'-dpng','-r300','-opengl');
        name = [outpath '/',RunID,'_phswt_',num2str(floor(step/nop))]; % figure 3
        print(fh3,name,'-dpng','-r300','-opengl');
        name = [outpath '/',RunID,'_vseg',num2str(floor(step/nop))]; % figure 4
        print(fh4,name,'-dpng','-r300','-opengl');
        name = [outpath '/',RunID,'_xc',num2str(floor(step/nop))]; % figure 5
        print(fh5,name,'-dpng','-r300','-opengl');
        name = [outpath '/',RunID,'_rho',num2str(floor(step/nop))]; % figure 6
        print(fh6,name,'-dpng','-r300','-opengl');
        name = [outpath '/',RunID,'_chm',num2str(floor(step/nop))]; % figure 7
        print(fh7,name,'-dpng','-r300','-opengl');
%         name = [outpath '/',RunID,'_xcd',num2str(floor(step/nop))]; % figure 8
%         print(fh8,name,'-dpng','-r300','-opengl');

    end
    
    name = [outpath '/',RunID,'_phsdg',num2str(floor(step/nop))]; % figure 9 phase diagram
    print(fh9,name,'-dpng','-r300','-opengl');
    name = [outpath '/',RunID,'_consv',num2str(floor(step/nop))]; % figure 10 conserved quantities
    print(fh10,name,'-dpng','-r300','-opengl');
    
    name = [outpath '/',RunID,'_',num2str(step/nop)];
    save(name,'U','W','P','Pt','xFe','xSi','cFe','cSi','csFe','clFe','csSi','clSi',...
        'fsFe','flFe','fsSi','flSi','S','XFe','XSi','CFe','CSi','FsFe','FlFe','FsSi','FlSi',...
        'phisFe','philFe','phisSi','philSi','rho','Eta','segsFe','seglFe','segsSi','seglSi',...
        'Ksgr_x','Ksgr_f','Ksgr_m','xP','zP','xU','zU');
    name = [outpath '/',RunID,'_cont'];
    save(name,'U','W','P','Pt','xFe','xSi','cFe','cSi','csFe','clFe','csSi','clSi',...
        'fsFe','flFe','fsSi','flSi','S','XFe','XSi','CFe','CSi','FsFe','FlFe','FsSi','FlSi',...
        'phisFe','philFe','phisSi','philSi','rho','Eta','segsFe','seglFe','segsSi','seglSi',...
        'Ksgr_x','Ksgr_f','Ksgr_m','xP','zP','xU','zU');
    name = [outpath,'/',RunID,'_hist'];
    save(name,'HIST');
    
    if step == 0
        logfile = [outpath '/',RunID,'.log'];
        if exist(logfile,'file'); delete(logfile); end
        diary(logfile)
    end
end

