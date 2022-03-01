% planetesimal: generate output

% print output header
fprintf(1,'\n*****  prepare output frame %d  *****\n',NUM.step/RUN.nop);


%% plot output figures
if RUN.plot
    
    % prepare for plotting
    TX = {'Interpreter','Latex'}; FS = {'FontSize',12};
    TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',10};
    UN = {'Units','Centimeters'};
    
    if NUM.Nx <= 10 && NUM.Nz <= 10  % create 0D plots
        
        fh1 = figure(1); clf;
        subplot(4,1,1)
        plot(HST.time/NUM.yr,HST.T(:,2),'LineWidth',2); axis xy tight; box on;
        title('$T [^\circ$C]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(4,1,2)
        plot(HST.time/NUM.yr,HST.xFe(:,2).*100,'LineWidth',2); axis xy tight; box on;
        title('$x_{Fe}$ [wt\% Fe-FeS]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(4,1,3)
        plot(HST.time/NUM.yr,HST.cFe(:,2).*100,'LineWidth',2); axis xy tight; box on;
        title('$\bar{c}_{Fe}$ [wt\% S]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(4,1,4)
        plot(HST.time/NUM.yr,HST.cSi(:,2).*100,'LineWidth',2); axis xy tight; box on;
        title('$\bar{c}_{Si}$ [wt\% SiO$_2$]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        xlabel('Time [yr]',TX{:},FS{:});
        
        fh2 = figure(2); clf;
        subplot(4,1,1)
        plot(HST.time/NUM.yr,HST.fFel(:,2).*100,'LineWidth',2); axis xy tight; box on;
        title(['$f_{Fe}$ [wt\%]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(4,1,2)
        plot(HST.time/NUM.yr,HST.fSil(:,2).*100,'LineWidth',2); axis xy tight; box on;
        title(['$f_{Si}$ [wt\%]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(4,1,3)
        plot(HST.time/NUM.yr,HST.rho(:,2),'LineWidth',2); axis xy tight; box on;
        title('$\bar{\rho}$ [kg/m$^3$]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(4,1,4)
        plot(HST.time/NUM.yr,log10(HST.eta(:,2)),'LineWidth',2); axis xy tight; box on;
        title('$\bar{\eta}$ [log$_{10}$ Pas]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        xlabel('Time [hr]',TX{:},FS{:});
        
    elseif NUM.Nx <= 10
        
        fh1 = figure(1); clf;
        subplot(1,4,1)
        plot(mean(SOL.T(2:end-1,2:end-1),2),NUM.zP(2:end-1),'LineWidth',2); axis ij tight; box on;
        title('$T [^\circ$C]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        ylabel('Depth [m]',TX{:},FS{:});
        subplot(1,4,2)
        plot(mean(CHM.xFe(2:end-1,2:end-1),2),NUM.zP(2:end-1),'LineWidth',2); axis ij tight; box on;
        title('$x_{Fe}$ [wt\% Fe-FeS]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,4,3)
        plot(mean(CHM.cFe(2:end-1,2:end-1),2),NUM.zP(2:end-1),'LineWidth',2); axis ij tight; box on;
        title('$\bar{c}_{Fe}$ [wt\% S]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,4,4)
        plot(mean(CHM.cSi(2:end-1,2:end-1),2),NUM.zP(2:end-1),'LineWidth',2); axis ij tight; box on;
        title('$\bar{c}_{Si}$ [wt\% SiO$_2$]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        
        fh2 = figure(2); clf;
        subplot(1,4,1)
        plot(mean(CHM.fFel(2:end-1,2:end-1),2),NUM.zP(2:end-1),'LineWidth',2); axis ij tight; box on;
        title(['$f_{Fe}^\ell$ [wt\%]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
        ylabel('Depth [m]',TX{:},FS{:});
        subplot(1,4,2)
        plot(mean(CHM.fSil(2:end-1,2:end-1),2),NUM.zP(2:end-1),'LineWidth',2); axis ij tight; box on;
        title(['$f_{Si}^\ell$ [wt\%]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,4,3)
        plot(mean(CHM.GFe(2:end-1,2:end-1)./MAT.rho(2:end-1,2:end-1).*100.*NUM.yr,2),NUM.zP(2:end-1),'LineWidth',2); axis ij tight; box on;
        title('$\Gamma_{Fe}$ [wt\%/yr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,4,4)
        plot(mean(CHM.GSi(2:end-1,2:end-1)./MAT.rho(2:end-1,2:end-1).*100.*NUM.yr,2),NUM.zP(2:end-1),'LineWidth',2); axis ij tight; box on;
        title('$\Gamma_{Si}$ [wt\%/yr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        
        fh3 = figure(3); clf;
        subplot(1,4,1)
        plot(mean(-SOL.W(:,2:end-1),2)*NUM.yr,NUM.zW,'LineWidth',2); axis ij tight; box on;
        title(['$W$ [m/yr]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
        ylabel('Depth [m]',TX{:},FS{:});
        subplot(1,4,2)
        plot(mean(SOL.P(2:end-1,2:end-1),2),NUM.zP(2:end-1),'LineWidth',2); axis ij tight; box on;
        title(['$P$ [Pa]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,4,3)
        plot(mean(MAT.rho(2:end-1,2:end-1),2),NUM.zP(2:end-1),'LineWidth',2); axis ij tight; box on;
        title('$\bar{\rho}$ [kg/m$^3$]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,4,4)
        plot(mean(log10(MAT.Eta(2:end-1,2:end-1)),2),NUM.zP(2:end-1),'LineWidth',2); axis ij tight; box on;
        title('$\bar{\eta}$ [log$_{10}$ Pas]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        
        fh4 = figure(4); clf;
        subplot(1,4,1)
        plot(mean(-(MAT.phiFel(1:end-1,2:end-1)+MAT.phiFel(2:end,2:end-1))/2.*segFel(:,2:end-1),2)*NUM.yr,NUM.zW,'LineWidth',2); axis ij tight; box on;
        title(['$w_{\Delta,Fe}^\ell$ [m/yr]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
        ylabel('Depth [m]',TX{:},FS{:});
        subplot(1,4,2)
        plot(mean(-(MAT.phiFes(1:end-1,2:end-1)+MAT.phiFes(2:end,2:end-1))/2.*segFes(:,2:end-1),2)*NUM.yr,NUM.zW,'LineWidth',2); axis ij tight; box on;
        title(['$w_{\Delta,Fe}^s$ [m/yr]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,4,3)
        plot(mean(-(MAT.phiSil(1:end-1,2:end-1)+MAT.phiSil(2:end,2:end-1))/2.*segSil(:,2:end-1),2)*NUM.yr,NUM.zW,'LineWidth',2); axis ij tight; box on;
        title(['$w_{\Delta,Si}^\ell$ [m/yr]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,4,4)
        plot(mean(-(MAT.phiSis(1:end-1,2:end-1)+MAT.phiSis(2:end,2:end-1))/2.*segSis(:,2:end-1),2)*NUM.yr,NUM.zW,'LineWidth',2); axis ij tight; box on;
        title(['$w_{\Delta,Si}^s$ [m/yr]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
        
    else
        
        % prepare for plotting
        TX = {'Interpreter','Latex'}; FS = {'FontSize',12};
        TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',10};
        UN = {'Units','Centimeters'};
        xq = round(NUM.N/10/2)+1:round(NUM.N/10):NUM.nxP;  % x-indexes for quiver plots
        zq = round(NUM.N/10/2)+1:round(NUM.N/10):NUM.nxP;  % z-indexes for quiver plots
        
        
        %% plot velocity and temperature solutions
        fh1  = figure(1); clf
        figure(1)
        subplot(2,3,1)
        imagesc(NUM.xU(:),NUM.zU(2:end-1),SOL.U(2:end-1,:)*NUM.yr); hold on;
        quiver(NUM.xP(xq),NUM.zP(zq),SOL.UP(zq,xq),SOL.WP(zq,xq),'k')
        colormap(subplot(2,3,1),cm2)
        axis ij equal tight;
        colorbar
        title('v^*_x [m/yr]')
        
        subplot(2,3,2)
        imagesc(NUM.xW(2:end-1),NUM.zW(:),-SOL.W(:,2:end-1)*NUM.yr); hold on;
        quiver(NUM.xP(xq),NUM.zP(zq),SOL.UP(zq,xq),SOL.WP(zq,xq),'k')
        colormap(subplot(2,3,2),cm2)
        axis ij equal tight;
        colorbar
        title('v_z^*-velocity [m/yr]')
        
        subplot(2,3,3)
        imagesc(NUM.xP(2:end-1),NUM.zP(2:end-1),SOL.P(2:end-1,2:end-1)); hold on;
        colormap(subplot(2,3,3),cm2)
        axis ij equal tight;
        colorbar
        title('dynamic pressure [Pa]')
        
        subplot(2,3,4)
        imagesc(NUM.xP(2:end-1),NUM.zP(2:end-1),SOL.T(2:end-1,2:end-1));
        colormap(subplot(2,3,4),flipud(cm1))
        axis ij equal tight;
        colorbar
        title('Temperature [C]')
        
        subplot(2,3,5);
        imagesc(NUM.xP(2:end-1),NUM.zP(2:end-1),MAT.rho(2:end-1,2:end-1));
        colormap(subplot(2,3,5),cm1)
        axis ij equal tight;
        colorbar
        title('bulk density [kgm^-^3]')
        
        subplot(2,3,6);
        imagesc(NUM.xP(2:end-1),NUM.zP(2:end-1),log10(MAT.Eta(2:end-1,2:end-1)));
        colormap(subplot(2,3,6),cm1)
        axis ij equal tight;
        colorbar
        title('\eta [log_1_0 Pas]')
        
        
        %% plot phase vol fractions
        fh2 = figure(2); clf;
        figure(2) % plot phase fractions
        subplot(2,2,1);
        imagesc(NUM.xP(2:end-1),NUM.zP(2:end-1),MAT.phiSil(2:end-1,2:end-1));
        colorbar
        axis ij equal tight;
        title('\phi_{Si}^l')
        subplot(2,2,2);
        imagesc(NUM.xP(2:end-1),NUM.zP(2:end-1),MAT.phiSis(2:end-1,2:end-1));
        colorbar
        axis ij equal tight;
        title('\phi_{Si}^s')
        subplot(2,2,3);
        imagesc(NUM.xP(2:end-1),NUM.zP(2:end-1),MAT.phiFel(2:end-1,2:end-1));
        colorbar
        axis ij equal tight;
        title('\phi_{Fe}^l')
        subplot(2,2,4);
        imagesc(NUM.xP(2:end-1),NUM.zP(2:end-1),MAT.phiFes(2:end-1,2:end-1));
        colorbar
        axis ij equal tight;
        title('\phi_{Fe}^s')
        
        
        %% plot phase
        fh3 = figure(3); clf
        figure(3) % plot phase fractions
        subplot(2,2,1);
        imagesc(NUM.xP(2:end-1),NUM.zP(2:end-1),CHM.fSil(2:end-1,2:end-1));
        colorbar
        axis ij equal tight
        title('f_{Si}^l [wt]')
        subplot(2,2,2);
        imagesc(NUM.xP(2:end-1),NUM.zP(2:end-1),CHM.fSis(2:end-1,2:end-1));
        colorbar
        axis ij equal tight
        title('f_{Si}^s [wt]')
        subplot(2,2,3);
        imagesc(NUM.xP(2:end-1),NUM.zP(2:end-1),CHM.fFel(2:end-1,2:end-1));
        colorbar
        axis ij equal tight
        title('f_{Fe}^l [wt]')
        subplot(2,2,4);
        imagesc(NUM.xP(2:end-1),NUM.zP(2:end-1),CHM.fFes(2:end-1,2:end-1));
        colorbar
        axis ij equal tight
        title('f_{Fe}^s [wt]')
        
        
        %% plot phase segregation
        fh4 = figure(4); clf
        figure(4) % plot phase segregation
        subplot(2,2,1);
        imagesc(NUM.xW(2:end-1),NUM.zW(:),-segSil(:,2:end-1)*NUM.yr);
        colorbar
        axis ij equal tight
        title('\Delta v_{Si}^l [m/yr]')
        
        subplot(2,2,2);
        imagesc(NUM.xW(2:end-1),NUM.zW(:),-segSis(:,2:end-1)*NUM.yr);
        colorbar
        axis ij equal tight
        title('\Delta v_{Si}^s [m/yr]')
        
        subplot(2,2,3);
        imagesc(NUM.xW(2:end-1),NUM.zW(:),-segFel(:,2:end-1)*NUM.yr);
        colorbar
        axis ij equal tight
        title('\Delta v_{Fe}^l [m/yr]')
        
        subplot(2,2,4);
        imagesc(NUM.xW(2:end-1),NUM.zW(:),-segFes(:,2:end-1)*NUM.yr);
        colorbar
        axis ij equal tight
        title('\Delta v_{Fe}^s [m/yr]')
        
        
        %% plot system fractions and chemical solutions
        fh5 = figure(5); clf
        figure(5)
        subplot(2,2,1)
        imagesc(NUM.xP(2:end-1),NUM.zP(2:end-1),CHM.xFe(2:end-1,2:end-1))
        colorbar
        axis ij equal tight
        title('x_{Fe} [wt]')
        
        subplot(2,2,2)
        imagesc(NUM.xP(2:end-1),NUM.zP(2:end-1),CHM.xSi(2:end-1,2:end-1))
        colorbar
        axis ij equal tight
        title('x_{Si} [wt]')
        
        subplot(2,2,3)
        imagesc(NUM.xP(2:end-1),NUM.zP(2:end-1),CHM.cFe(2:end-1,2:end-1))
        colorbar
        axis ij equal tight
        title('c_{Fe} [wt]')
        
        subplot(2,2,4)
        imagesc(NUM.xP(2:end-1),NUM.zP(2:end-1),CHM.cSi(2:end-1,2:end-1))
        colorbar
        axis ij equal tight
        title('c_{Si} [wt]')
        
        
        %% plot phase densities
        fh6 = figure(6); clf
        figure(6)
        subplot(2,2,1)
        imagesc(NUM.xP(2:end-1),NUM.zP(2:end-1),MAT.rhoSil(2:end-1,2:end-1))
        colorbar
        axis ij equal tight
        title('\rho_{Si}^l')
        
        subplot(2,2,2)
        imagesc(NUM.xP(2:end-1),NUM.zP(2:end-1),MAT.rhoSis(2:end-1,2:end-1))
        colorbar
        axis ij equal tight
        title('\rho_{Si}^s')
        
        subplot(2,2,3)
        imagesc(NUM.xP(2:end-1),NUM.zP(2:end-1),MAT.rhoFel(2:end-1,2:end-1))
        colorbar
        axis ij equal tight
        title('\rho_{Fe}^l')
        
        subplot(2,2,4)
        imagesc(NUM.xP(2:end-1),NUM.zP(2:end-1),MAT.rhoFes(2:end-1,2:end-1))
        colorbar
        axis ij equal tight
        title('\rho_{Fe}^s')
        
        
        %% plot chemical phase compositions
        fh7 = figure(7); clf
        figure(7)
        subplot(2,2,1)
        imagesc(NUM.xP(2:end-1),NUM.zP(2:end-1),CHM.clSi(2:end-1,2:end-1))
        colorbar
        axis ij equal tight
        title('c_{Si}^l')
        
        subplot(2,2,2)
        imagesc(NUM.xP(2:end-1),NUM.zP(2:end-1),CHM.csSi(2:end-1,2:end-1))
        colorbar
        axis ij equal tight
        title('c_{Si}^s')
        
        subplot(2,2,3)
        imagesc(NUM.xP(2:end-1),NUM.zP(2:end-1),CHM.clFe(2:end-1,2:end-1))
        colorbar
        axis ij equal tight
        title('c_{Fe}^l')
        
        subplot(2,2,4)
        imagesc(NUM.xP(2:end-1),NUM.zP(2:end-1),CHM.csFe(2:end-1,2:end-1))
        colorbar
        axis ij equal tight
        title('c_{Fe}^s')
        
        
        %% plot Partial densities and
        fh8 = figure(8); clf;
        figure(8)
        subplot(2,2,1)
        imagesc(NUM.xP(2:end-1),NUM.zP(2:end-1),CHM.CFe(2:end-1,2:end-1))
        colorbar
        axis ij equal tight
        title('C_{Fe} ')
        
        subplot(2,2,2)
        imagesc(NUM.xP(2:end-1),NUM.zP(2:end-1),CHM.CSi(2:end-1,2:end-1))
        colorbar
        axis ij equal tight
        title('C_{Si} ')
        
        subplot(2,2,3)
        imagesc(NUM.xP(2:end-1),NUM.zP(2:end-1),CHM.XFe(2:end-1,2:end-1))
        colorbar
        axis ij equal tight
        title('X_{Fe}')
        
        subplot(2,2,4)
        imagesc(NUM.xP(2:end-1),NUM.zP(2:end-1),CHM.XSi(2:end-1,2:end-1))
        colorbar
        axis ij equal tight
        title('X_{Si}')
        
    end
    
    %% plot phase diagrams
    fh9 = figure(9); clf;
    TT = linspace(CHM.TSi1,CHM.TSi2,1e3);
    cc = [linspace(CHM.cphsSi2,(CHM.perCsSi+CHM.perClSi)/2,round((CHM.perTSi-CHM.TSi1)./(CHM.TSi2-CHM.TSi1)*1e3)),linspace((CHM.perCsSi+CHM.perClSi)/2,CHM.cphsSi1,round((CHM.perTSi-CHM.TSi2)./(CHM.TSi1-CHM.TSi2)*1e3))];
    [~,CCxSi,CClSi] = equilibrium(TT,cc,0.*TT,CHM.TSi1,CHM.TSi2,CHM.cphsSi1,CHM.cphsSi2,...
                                  CHM.perTSi,CHM.perCsSi,CHM.perClSi,CHM.clap,CHM.PhDgSi,TINY);
    
    TT2 = linspace(CHM.TFe1,CHM.TFe2,1000);
    cc2 = linspace(CHM.cphsFe2,CHM.cphsFe1,length(TT2));
    [~,CCxFe,CClFe] = equilibrium(TT2,cc2,0.*TT2,CHM.TFe1,CHM.TFe2,CHM.cphsFe1,CHM.cphsFe2,...
                                  CHM.perTFe,CHM.perCsFe,CHM.perClFe,CHM.clap,CHM.PhDgFe,TINY);
    subplot(1,2,1)
    plot(CCxSi,TT,'k-','LineWidth',2); axis tight; hold on; box on;
    plot(CClSi,TT,'k-','LineWidth',2); axis tight; hold on; box on;
    plot(CHM.cSi(2:end-1),SOL.T(2:end-1), '.k', CHM.csSi(2:end-1),SOL.T(2:end-1), '.b', CHM.clSi(2:end-1),SOL.T(2:end-1), '.r','MarkerSize',15)
    
    subplot(1,2,2)
    plot(CCxFe,TT2,'k-','LineWidth',2); axis tight; hold on; box on;
    plot(CClFe,TT2,'k-','LineWidth',2); axis tight; hold on; box on;
    plot(CHM.cFe(2:end-1),SOL.T(2:end-1), '.k', CHM.csFe(2:end-1),SOL.T(2:end-1), '.b', CHM.clFe(2:end-1),SOL.T(2:end-1), '.r','MarkerSize',15)
    
    drawnow;
    
    
    %% plot conserved quantities
    fh10 = figure(10); clf;
    % mass/energy conservation
    subplot(5,1,1)
    plot(HST.time./NUM.yr,HST.EM  ,'k-','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
    ylabel('consv. $Mass$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(5,1,2)
    plot(HST.time./NUM.yr,HST.EH  ,'k-','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
    ylabel('consv. $H$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(5,1,3)
    plot(HST.time./NUM.yr,HST.EXFe,'k-','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
    ylabel('consv. $X_{Fe}$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(5,1,4)
    plot(HST.time./NUM.yr,HST.ECFe,'k-','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
    ylabel('consv. $C_{Fe}$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(5,1,5)
    plot(HST.time./NUM.yr,HST.ECSi,'k-','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
    ylabel('consv. $C_{Si}$',TX{:},FS{:}); set(gca,TL{:},TS{:})
    xlabel('Time [yr]',TX{:},FS{:});
    
end


%% save output
if RUN.save
    if NUM.Nx <= 10 && NUM.Nz <= 10  % save 0D plots
        
        name = [outpath '/',RUN.ID,'_txc_',num2str(floor(NUM.step/RUN.nop))]; % figure 1
        print(fh1,name,'-dpng','-r300','-opengl');
        name = [outpath '/',RUN.ID,'_fre_',num2str(floor(NUM.step/RUN.nop))]; % figure 2
        print(fh2,name,'-dpng','-r300','-opengl');
        
    elseif NUM.Nx <= 10   % save 1D plots
        
        name = [outpath '/',RUN.ID,'_txc_',num2str(floor(NUM.step/RUN.nop))]; % figure 1
        print(fh1,name,'-dpng','-r300','-opengl');
        name = [outpath '/',RUN.ID,'_fmr_',num2str(floor(NUM.step/RUN.nop))]; % figure 2
        print(fh2,name,'-dpng','-r300','-opengl');
        name = [outpath '/',RUN.ID,'_vp_',num2str(floor(NUM.step/RUN.nop))]; % figure 2
        print(fh2,name,'-dpng','-r300','-opengl');
        name = [outpath '/',RUN.ID,'_seg_',num2str(floor(NUM.step/RUN.nop))]; % figure 2
        print(fh2,name,'-dpng','-r300','-opengl');
        
    else  % save 2D plots
        
        name = [outpath '/',RUN.ID,'_vp_',num2str(floor(NUM.step/RUN.nop))]; % figure 1
        print(fh1,name,'-dpng','-r300','-opengl');
        name = [outpath '/',RUN.ID,'_phsvol_',num2str(floor(NUM.step/RUN.nop))]; % figure 2
        print(fh2,name,'-dpng','-r300','-opengl');
        name = [outpath '/',RUN.ID,'_phswt_',num2str(floor(NUM.step/RUN.nop))]; % figure 3
        print(fh3,name,'-dpng','-r300','-opengl');
        name = [outpath '/',RUN.ID,'_vseg',num2str(floor(NUM.step/RUN.nop))]; % figure 4
        print(fh4,name,'-dpng','-r300','-opengl');
        name = [outpath '/',RUN.ID,'_xcb',num2str(floor(NUM.step/RUN.nop))]; % figure 5
        print(fh5,name,'-dpng','-r300','-opengl');
        name = [outpath '/',RUN.ID,'_rho',num2str(floor(NUM.step/RUN.nop))]; % figure 5
        print(fh6,name,'-dpng','-r300','-opengl');
        name = [outpath '/',RUN.ID,'_chm',num2str(floor(NUM.step/RUN.nop))]; % figure 5
        print(fh7,name,'-dpng','-r300','-opengl');
        name = [outpath '/',RUN.ID,'_xcd',num2str(floor(NUM.step/RUN.nop))]; % figure 5
        print(fh8,name,'-dpng','-r300','-opengl');

    end
    
    name = [outpath '/',RUN.ID,'_phsdg',num2str(floor(NUM.step/RUN.nop))]; % figure 9 phase diagram
    print(fh9,name,'-dpng','-r300','-opengl');
    name = [outpath '/',RUN.ID,'_consv',num2str(floor(NUM.step/RUN.nop))]; % figure 10 conserved quantities
    print(fh10,name,'-dpng','-r300','-opengl');
    
    name = [outpath '/',RUN.ID,'_',num2str(NUM.step/RUN.nop)];
    save(name,'NUM','PHY','SOL','MAT','CHM','HST');
    name = [outpath '/',RUN.ID,'_cont'];
    save(name,'NUM','PHY','SOL','MAT','CHM','HST');
    
    if NUM.step == 0
        logfile = [outpath '/',RUN.ID,'.log'];
        if exist(logfile,'file'); delete(logfile); end
        diary(logfile)
    end
end

