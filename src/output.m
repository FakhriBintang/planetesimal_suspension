% planetesimal: generate output

% print output header
fprintf(1,'\n*****  prepare output frame %d  *****\n',step/nop);


%% plot output figures
if plot_op

    % prepare for plotting
    % prepare for plotting
    TX = {'Interpreter','Latex'}; FS = {'FontSize',12};
    TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',10};
    UN = {'Units','Centimeters'};
    CL = {'Color',[0.0 0.0 0.0],[0.80 0.15 0.10],[0.10 0.15 0.65],[0.45 0.60 0.95]};
    LW = {'LineWidth',2};

    if Nx <= 10 && Nz <= 10  % create 0D plots

        fh1 = figure(1); clf;
        subplot(4,1,1)
        plot(HST.time/yr,HST.T,'LineWidth',2); axis xy tight; box on;
        title('$T [^\circ$k]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(4,1,2)
        plot(HST.time/yr,HST.xFe.*100,'LineWidth',2); axis xy tight; box on;
        title('$x_{Fe}$ [wt\% Fe-FeS]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(4,1,3)
        plot(HST.time/yr,HST.cFe.*100,'LineWidth',2); axis xy tight; box on;
        title('$\bar{c}_{Fe}$ [wt\% S]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(4,1,4)
        plot(HST.time/yr,HST.cSi.*100,'LineWidth',2); axis xy tight; box on;
        title('$\bar{c}_{Si}$ [wt\% SiO$_2$]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        xlabel('Time [yr]',TX{:},FS{:});

        fh2 = figure(2); clf;
        subplot(4,1,1)
        plot(HST.time/yr,HST.flFe.*100,'LineWidth',2); axis xy tight; box on;
        title('$f_{Fe}$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(4,1,2)
        plot(HST.time/yr,HST.flSi.*100,'LineWidth',2); axis xy tight; box on;
        title('$f_{Si}$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(4,1,3)
        plot(HST.time/yr,HST.rho,'LineWidth',2); axis xy tight; box on;
        title('$\bar{\rho}$ [kg/m$^3$]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(4,1,4)
        plot(HST.time/yr,log10(HST.Eta),'LineWidth',2); axis xy tight; box on;
        title('$\bar{\eta}$ [log$_{10}$ Pas]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        xlabel('Time [yr]',TX{:},FS{:});

    elseif Nx <= 10 % 1D plots

        fh1 = figure(1); clf;
        subplot(1,4,1)
        plot(mean(T(2:end-1,2:end-1),2),zP(2:end-1)./1000,'LineWidth',2); axis ij tight; box on;
        title('$T [^\circ$K]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        ylabel('Depth [km]',TX{:},FS{:});
        subplot(1,4,2)
        plot(mean(xFe(2:end-1,2:end-1),2).*100,zP(2:end-1)./1000,'LineWidth',2); hold on;  axis ij tight; box on;
        plot(mean(xSi(2:end-1,2:end-1),2).*100,zP(2:end-1)./1000,'LineWidth',2);
        title('$x_{Fe}$ / $x_{Si}$ [wt\% Fe-FeS]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        legend('$x_{\mathrm{Fe}}$', '$x_{\mathrm{Si}}$',TX{:})
        subplot(1,4,3)
        plot(mean(clFe(2:end-1,2:end-1),2).*100,zP(2:end-1)./1000,'-r','LineWidth',2); axis ij tight; box on; hold on
        plot(mean(csFe(2:end-1,2:end-1),2).*100,zP(2:end-1)./1000,'-b','LineWidth',2);
        plot(mean(cFe(2:end-1,2:end-1),2).*100,zP(2:end-1)./1000,'k-','LineWidth',2);
        title('$\bar{c}_{Fe}$ [wt\% S]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,4,4)
        plot((mean(clSi(2:end-1,2:end-1),2)+cSimin).*100,zP(2:end-1)./1000,'-r','LineWidth',2); axis ij tight; box on; hold on
        plot((mean(csSi(2:end-1,2:end-1),2)+cSimin).*100,zP(2:end-1)./1000,'-b','LineWidth',2);
        plot((mean(cSi(2:end-1,2:end-1),2)+cSimin).*100,zP(2:end-1)./1000,'k-','LineWidth',2);
        title('$\bar{c}_{Si}$ [wt\% SiO$_2$]',TX{:},FS{:}); set(gca,TL{:},TS{:});

        fh2 = figure(2); clf;
        subplot(1,4,1)
        plot(mean(flFe(2:end-1,2:end-1),2).*100,zP(2:end-1)./1000,'LineWidth',2); axis ij tight; box on; hold on
        plot(mean(flSi(2:end-1,2:end-1),2).*100,zP(2:end-1)./1000,'LineWidth',2);
        title('$f_{j}^\ell$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        legend('Fe', 'Si')
        ylabel('Depth [km]',TX{:},FS{:});
        subplot(1,4,2)
        plot(mean(phisFe(2:end-1,2:end-1),2).*100,zP(2:end-1)./1000,'LineWidth',2); axis ij tight; box on; hold on
        plot(mean(philFe(2:end-1,2:end-1),2).*100,zP(2:end-1)./1000,'LineWidth',2);
        plot(mean(phisSi(2:end-1,2:end-1),2).*100,zP(2:end-1)./1000,'LineWidth',2);
        plot(mean(philSi(2:end-1,2:end-1),2).*100,zP(2:end-1)./1000,'LineWidth',2);
        title('$\phi_{j}^i$ [vol\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        legend('\phi_{Fe}^s','\phi_{Fe}^{l}','\phi_{Si}^s','\phi_{Si}^{l}')
        subplot(1,4,3)
        plot(mean(GFes(2:end-1,2:end-1)./rho(2:end-1,2:end-1).*100.*yr,2),zP(2:end-1)./1000,'LineWidth',2); axis ij tight; box on;
        title('$\Gamma_{Fe}^s$ [wt\%/yr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,4,4)
        plot(mean(GSis(2:end-1,2:end-1)./rho(2:end-1,2:end-1).*100.*yr,2),zP(2:end-1)./1000,'LineWidth',2); axis ij tight; box on;
        title('$\Gamma_{Si}^s$ [wt\%/yr]',TX{:},FS{:}); set(gca,TL{:},TS{:});

        fh3 = figure(3); clf;
        subplot(1,4,1)
        plot(mean(-(phisFe(1:end-1,2:end-1)+phisFe(2:end,2:end-1))/2.*segsFe(:,2:end-1),2)*yr,zW,'LineWidth',2); hold on
        plot(mean(-(philFe(1:end-1,2:end-1)+philFe(2:end,2:end-1))/2.*seglFe(:,2:end-1),2)*yr,zW,'LineWidth',2);
        plot(mean(-(phisSi(1:end-1,2:end-1)+phisSi(2:end,2:end-1))/2.*segsSi(:,2:end-1),2)*yr,zW,'LineWidth',2);
        plot(mean(-(philSi(1:end-1,2:end-1)+philSi(2:end,2:end-1))/2.*seglSi(:,2:end-1),2)*yr,zW,'LineWidth',2);
        plot(mean(-W(:,2:end-1),2)*yr,zW,'k','LineWidth',2); axis ij tight; box on;
        title('$W$ [m/yr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        ylabel('Depth [km]',TX{:},FS{:});
        subplot(1,4,2)
        plot(mean(P(2:end-1,2:end-1),2),zP(2:end-1)./1000,'LineWidth',2); axis ij tight; box on;
        title('$P$ [Pa]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,4,3)
        plot(mean(rho(2:end-1,2:end-1),2),zP(2:end-1)./1000,'LineWidth',2); axis ij tight; box on;
        title('$\bar{\rho}$ [kg/m$^3$]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,4,4)
        plot(mean(log10(Eta(2:end-1,2:end-1)),2),zP(2:end-1)./1000,'LineWidth',2); axis ij tight; box on;
        title('$\bar{\eta}$ [log$_{10}$ Pas]',TX{:},FS{:}); set(gca,TL{:},TS{:});

        fh4 = figure(4); clf;
        subplot(1,4,1)
        plot(mean(-(philFe(1:end-1,2:end-1)+philFe(2:end,2:end-1))/2.*seglFe(:,2:end-1),2)*yr,zW,'LineWidth',2); axis ij tight; box on;
        title('$w_{\Delta,Fe}^\ell$ [m/yr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        ylabel('Depth [km]',TX{:},FS{:});
        subplot(1,4,2)
        plot(mean(-(phisFe(1:end-1,2:end-1)+phisFe(2:end,2:end-1))/2.*segsFe(:,2:end-1),2)*yr,zW,'LineWidth',2); axis ij tight; box on;
        title('$w_{\Delta,Fe}^s$ [m/yr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,4,3)
        plot(mean(-(philSi(1:end-1,2:end-1)+philSi(2:end,2:end-1))/2.*seglSi(:,2:end-1),2)*yr,zW,'LineWidth',2); axis ij tight; box on;
        title('$w_{\Delta,Si}^\ell$ [m/yr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,4,4)
        plot(mean(-(phisSi(1:end-1,2:end-1)+phisSi(2:end,2:end-1))/2.*segsSi(:,2:end-1),2)*yr,zW,'LineWidth',2); axis ij tight; box on;
        title('$w_{\Delta,Si}^s$ [m/yr]',TX{:},FS{:}); set(gca,TL{:},TS{:});

        % plot conserved quantities
        fh5 = figure(5); clf;
        subplot(1,5,1)
        plot(mean(XFe(2:end-1,2:end-1),2),zP(2:end-1)./1000,'LineWidth',2); axis ij tight; box on; hold on;
        plot(mean(XSi(2:end-1,2:end-1),2),zP(2:end-1)./1000,'LineWidth',2);
        plot(mean(rho(2:end-1,2:end-1),2),zP(2:end-1)./1000,'-k','LineWidth',2);
        plot(mean(XFe(2:end-1,2:end-1),2)+ mean(XSi(2:end-1,2:end-1),2),zP(2:end-1)./1000,'--','LineWidth',2);
        legend('XFe', 'XSi', 'rho', 'XFe +XSi')
        ylabel('Depth [km]',TX{:},FS{:});
        subplot(1,5,2)
        plot(mean(CFe(2:end-1,2:end-1),2),zP(2:end-1)./1000,'LineWidth',2); axis ij tight; box on;
        title('$C_{Fe}$ ',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,5,3)
        plot(mean(CSi(2:end-1,2:end-1),2),zP(2:end-1)./1000,'LineWidth',2); axis ij tight; box on;
        title('$C_{Si}$ ',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,5,4)
        plot(mean(FlFe(2:end-1,2:end-1),2),zP(2:end-1)./1000,'LineWidth',2); axis ij tight; box on; hold on;
        plot(mean(FsFe(2:end-1,2:end-1),2),zP(2:end-1)./1000,'LineWidth',2);
        plot(mean(FlSi(2:end-1,2:end-1),2),zP(2:end-1)./1000,'LineWidth',2);
        plot(mean(FsSi(2:end-1,2:end-1),2),zP(2:end-1)./1000,'LineWidth',2);
        legend('FlFe','FsFe' , 'FlSi', 'FsSi')
        ylabel('Depth [km]',TX{:},FS{:});
        subplot(1,5,5)
        plot(mean(S(2:end-1,2:end-1),2),zP(2:end-1)./1000,'LineWidth',2); axis ij tight; box on;
        title('$S$ ',TX{:},FS{:}); set(gca,TL{:},TS{:});

    else
        %% prepare 2D plots
        % prepare for plotting
        TX = {'Interpreter','Latex'}; FS = {'FontSize',12};
        TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',10};
        UN = {'Units','Centimeters'};
        xq = round(N/10/2)+1:round(N/10):nxP;  % x-indexes for quiver plots
        zq = round(N/10/2)+1:round(N/10):nxP;  % z-indexes for quiver plots
        
        % set axis and border dimensions
        axh = 6.00*sqrt(D/L); axw = 6.00*sqrt(L/D)+1.50;
        ahs = 0.60; avs = 0.80;
        axb = 1.20; axt = 1.50;
        axl = 1.20; axr = 0.60;
        % plot velocity and temperature solutions
        if ~exist('fh1','var'); fh1 = figure(1);
        else; set(0, 'CurrentFigure', fh1); clf;
        end
        colormap(ocean);
        fh = axb + 2*axh + 1*avs + axt;
        fw = axl + 3*axw + 2*ahs + axr;
        set(fh1,UN{:},'Position',[13 13 fw fh]);
        set(fh1,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
        set(fh1,'Color','w','InvertHardcopy','off','Resize','off');
        ax(11) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
        ax(12) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
        ax(13) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+1*axh+1*avs axw axh]);
        ax(14) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
        ax(15) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
        ax(16) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);

        set(0,'CurrentFigure',fh1)
        set(fh1,'CurrentAxes',ax(11));
        imagesc(xU(:),zU(2:end-1),U(2:end-1,:)*yr); hold on; axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('$v^*_x$ [m/yr]',TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [km]',TX{:},FS{:});
        set(fh1,'CurrentAxes',ax(12));
        imagesc(xW(2:end-1),zW(:),-W(:,2:end-1)*yr); hold on; axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('$v^*_z$ [m/yr]',TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
        set(fh1,'CurrentAxes',ax(13));
        imagesc(xP(2:end-1),zP(2:end-1)./1000,P(2:end-1,2:end-1)); hold on; axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('dynamic pressure [Pa]',TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
        set(fh1,'CurrentAxes',ax(14));
        imagesc(xP(2:end-1),zP(2:end-1)./1000,T(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('T [$^{\circ}$k]',TX{:},FS{:}); ylabel('Depth [km]',TX{:},FS{:});
        set(fh1,'CurrentAxes',ax(15));
        imagesc(xP(2:end-1),zP(2:end-1)./1000,rho(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('$\bar{\rho} [\mathrm{kg m^{-3}}]$',TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel('Width [m]',TX{:},FS{:});
        set(fh1,'CurrentAxes',ax(16));
        imagesc(xP(2:end-1),zP(2:end-1)./1000,log10(Eta(2:end-1,2:end-1))); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('Log10($\bar{\eta})$ [Pas]',TX{:},FS{:}); set(gca,'YTickLabel',[]);
        sgtitle(['time = ',num2str(time/yr,3),' [yr]'],TX{:},FS{:},'Color','k');

        %% plot phase vol fractions
        if ~exist('fh2','var'); fh2 = figure(2);
        else; set(0, 'CurrentFigure', fh2); clf;
        end
        colormap(ocean);
        fh = axb + 2*axh + 1*avs + axt;
        fw = axl + 2*axw + 1*ahs + axr;
        set(fh2,UN{:},'Position',[1 1 fw fh]);
        set(fh2,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
        set(fh2,'Color','w','InvertHardcopy','off','Resize','off');
        ax(21) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
        ax(22) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
        ax(23) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
        ax(24) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);

        set(0,'CurrentFigure',fh2)
        set(fh2,'CurrentAxes',ax(21));
        imagesc(xP(2:end-1),zP(2:end-1)./1000,philSi(2:end-1,2:end-1)*100); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('$\phi_{Si}^l$ [vol\%]',TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [km]',TX{:},FS{:});
        set(fh2,'CurrentAxes',ax(22));
        imagesc(xP(2:end-1),zP(2:end-1)./1000,phisSi(2:end-1,2:end-1*100)); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('$\phi_{Si}^s$ [vol\%]',TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
        set(fh2,'CurrentAxes',ax(23));
        imagesc(xP(2:end-1),zP(2:end-1)./1000,philFe(2:end-1,2:end-1)*100); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('$\phi_{Fe}^l$ [vol\%]',TX{:},FS{:}); ylabel('Depth [km]',TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:});
        set(fh2,'CurrentAxes',ax(24));
        imagesc(xP(2:end-1),zP(2:end-1)./1000,phisFe(2:end-1,2:end-1)*100); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('$\phi_{Fe}^s$ [vol\%]',TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:}); set(gca,'YTickLabel',[]);
        sgtitle(['time = ',num2str(time/yr,3),' [yr]'],TX{:},FS{:},'Color','k');

        %% plot phase
if ~exist('fh3','var'); fh3 = figure(3);
    else; set(0, 'CurrentFigure', fh3); clf;
    end 
    colormap(ocean);
    fh = axb + 2*axh + 1*avs + axt;
    fw = axl + 2*axw + 1*ahs + axr;
    set(fh3,UN{:},'Position',[5 5 fw fh]);
    set(fh3,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh3,'Color','w','InvertHardcopy','off','Resize','off');
    ax(31) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
    ax(32) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
    ax(33) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(34) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
    
    set(0,'CurrentFigure',fh3)
    set(fh3,'CurrentAxes',ax(31));
    imagesc(xP(2:end-1),zP(2:end-1)./1000,flSi(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('$f_{Si}^l$ [wt\%]',TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [km]',TX{:},FS{:}); 
    set(fh3,'CurrentAxes',ax(32));
    imagesc(xP(2:end-1),zP(2:end-1)./1000,fsSi(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('$f_{Si}^s$ [wt\%]',TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh3,'CurrentAxes',ax(33));
     imagesc(xP(2:end-1),zP(2:end-1)./1000,flFe(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('$f_{Fe}^l$ [wt\%]',TX{:},FS{:}); ylabel('Depth [km]',TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:});
    set(fh3,'CurrentAxes',ax(34));
    imagesc(xP(2:end-1),zP(2:end-1)./1000,fsFe(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('$f_{Fe}^s$ [wt\%]',TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:}); set(gca,'YTickLabel',[]);
    sgtitle(['time = ',num2str(time/yr,3),' [yr]'],TX{:},FS{:},'Color','k');

    if ~exist('fh4','var'); fh4 = figure(4);
    else; set(0, 'CurrentFigure', fh4); clf;
    end 
    colormap(ocean);
    fh = axb + 2*axh + 1*avs + axt;
    fw = axl + 2*axw + 1*ahs + axr;
    set(fh4,UN{:},'Position',[7 7 fw fh]);
    set(fh4,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh4,'Color','w','InvertHardcopy','off','Resize','off');
    ax(41) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
    ax(42) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
    ax(43) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(44) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);

    set(0,'CurrentFigure',fh4)
    set(fh4,'CurrentAxes',ax(41));
    imagesc(xW(2:end-1),zW(:),-seglSi(:,2:end-1)*yr); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('$\Delta v_{Si}^l$ [m/yr]',TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [km]',TX{:},FS{:}); 
    set(fh4,'CurrentAxes',ax(42));
    imagesc(xW(2:end-1),zW(:),-segsSi(:,2:end-1)*yr); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('$\Delta v_{Si}^s$ [m/yr]',TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh4,'CurrentAxes',ax(43));
    imagesc(xW(2:end-1),zW(:),-seglFe(:,2:end-1)*yr); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('$\Delta v_{Fe}^l$ [m/yr]',TX{:},FS{:}); ylabel('Depth [km]',TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:});
    set(fh4,'CurrentAxes',ax(44));
    imagesc(xW(2:end-1),zW(:),-segsFe(:,2:end-1)*yr); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('$\Delta v_{Fe}^s$ [m/yr]',TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:}); set(gca,'YTickLabel',[]);
    sgtitle(['time = ',num2str(time/yr,3),' [yr]'],TX{:},FS{:},'Color','k');
    
    if ~exist('fh5','var'); fh5 = figure(5);
    else; set(0, 'CurrentFigure', fh5); clf;
    end 
    colormap(ocean);
    fh = axb + 2*axh + 1*avs + axt;
    fw = axl + 2*axw + 1*ahs + axr;
    set(fh5,UN{:},'Position',[9 9 fw fh]);
    set(fh5,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh5,'Color','w','InvertHardcopy','off','Resize','off');
    ax(51) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
    ax(52) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
    ax(53) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(54) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);

    set(0,'CurrentFigure',fh5)
    set(fh5,'CurrentAxes',ax(51));
    imagesc(xP(2:end-1),zP(2:end-1)./1000,xFe(2:end-1,2:end-1)*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('$x_{Fe}$ [wt\%]',TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [km]',TX{:},FS{:});
    set(fh5,'CurrentAxes',ax(52));
     imagesc(xP(2:end-1),zP(2:end-1)./1000,xSi(2:end-1,2:end-1)*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('$x_{Si}$ [wt\%]',TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh5,'CurrentAxes',ax(53));
     imagesc(xP(2:end-1),zP(2:end-1)./1000,cFe(2:end-1,2:end-1)*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('$c_{Fe}$ [wt\%]',TX{:},FS{:}); ylabel('Depth [km]',TX{:},FS{:});
    set(fh5,'CurrentAxes',ax(54));
    imagesc(xP(2:end-1),zP(2:end-1)./1000,(cSi(2:end-1,2:end-1)+cSimin)*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('$c_{Si}$ [wt\%]',TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel('Width [m]',TX{:},FS{:});
    sgtitle(['time = ',num2str(time/yr,3),' [yr]'],TX{:},FS{:},'Color','k');

    end

    %% plot phase diagrams
    fh9 = figure(9); clf;
    TT = linspace(TSi1,TSi2,1e3);
    cc = [linspace(cphsSi2,(percsSi+perclSi)/2,round((perTSi-TSi1)./(TSi2-TSi1)*1e3)),linspace((percsSi+perclSi)/2,cphsSi1,round((perTSi-TSi2)./(TSi1-TSi2)*1e3))];
    [~,CCxSi,CClSi] = equilibrium(TT,cc,0.*TT,TSi1,TSi2,cphsSi1,cphsSi2,...
        perTSi,percsSi,perclSi,clap,PhDgSi);

    TT2 = linspace(TFe1,TFe2,1000);
    cc2 = linspace(cphsFe2,cphsFe1,length(TT2));
    [~,CCxFe,CClFe] = equilibrium(TT2,cc2,0.*TT2,TFe1,TFe2,cphsFe1,cphsFe2,...
        perTFe,percsFe,perclFe,clap,PhDgFe);
    subplot(1,2,1)
    plot(CCxSi+cSimin,TT,'k-','LineWidth',2); axis tight; axis square; hold on; box on;
    plot(CClSi+cSimin,TT,'k-','LineWidth',2); axis tight; hold on; axis square; box on;
    plot(cSi(2:end-1,2:end-1)+cSimin,T(2:end-1,2:end-1), '.k', csSi(2:end-1,2:end-1)+cSimin,T(2:end-1,2:end-1), '.b', clSi(2:end-1,2:end-1)+cSimin,T(2:end-1,2:end-1), '.r','MarkerSize',25)
    ylim([TSi1 TSi2])
    title('SiO$_2$ Phase Diagram','Interpreter','latex','FontSize',18)
    xlabel('Major component [wt\% SiO$_2$]','Interpreter','latex','FontSize',15)
    ylabel('Temperature [$^\circ$k]','Interpreter','latex','FontSize',15)

    subplot(1,2,2)
    plot(CCxFe,TT2,'k-','LineWidth',2); axis tight; axis square; hold on; box on;
    plot(CClFe,TT2,'k-','LineWidth',2); axis tight; axis square; hold on; box on;
    plot(cFe(2:end-1,2:end-1),T(2:end-1,2:end-1), '.k', csFe(2:end-1,2:end-1),T(2:end-1,2:end-1), '.b', clFe(2:end-1,2:end-1),T(2:end-1,2:end-1), '.r','MarkerSize',25)
    ylim([TFe1 TSi2])
    title('Fe-FeS Phase Diagram','Interpreter','latex','FontSize',18)
    xlabel('Major component [wt\% S]','Interpreter','latex','FontSize',15)
    ylabel('Temperature [$^\circ$k]','Interpreter','latex','FontSize',15)
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

    %% temporary S diagnosis
    fh11 = figure(11); clf;
    Ssum = FlFe.*slFe + FsFe.*ssFe + FlSi.*slSi + FsSi.*ssSi; 
    subplot(1,6,1)
    plot(mean(S(2:end-1,2:end-1),2),zP(2:end-1)./1000,'LineWidth',2); axis ij tight; box on; hold on;
    plot(mean(Ssum(2:end-1,2:end-1),2),zP(2:end-1)./1000,'--k','LineWidth',2);
    title('S')
    subplot(1,6,2)
    plot(dSdt,zP(2:end-1)./1000,'LineWidth',2); axis ij tight; box on;
    title('dSdt')   
    subplot(1,6,3)
    plot(advn_S,zP(2:end-1)./1000,'LineWidth',2); axis ij tight; box on;
    title('advection S')
    subplot(1,6,4)
    plot(diss_T,zP(2:end-1)./1000,'LineWidth',2); axis ij tight; box on;
    title('diss T')
    subplot(1,6,5)
    plot(mean(slFe(2:end-1,2:end-1),2),zP(2:end-1)./1000,'LineWidth',2); axis ij tight; box on; hold on
    plot(mean(ssFe(2:end-1,2:end-1),2),zP(2:end-1)./1000,'LineWidth',2);
    plot(mean(slSi(2:end-1,2:end-1),2),zP(2:end-1)./1000,'LineWidth',2);
    plot(mean(ssSi(2:end-1,2:end-1),2),zP(2:end-1)./1000,'LineWidth',2);
    legend('s_{Fe}^l','s_{Fe}^s','s_{Si}^l','s_{Si}^s')
    title('phase entropy')  
    subplot(1,6,6)
    plot(mean(T(2:end-1,2:end-1),2),zP(2:end-1)./1000,'LineWidth',2); axis ij tight; box on; hold on;
    title('T')
end


%% save output
if save_op
    if Nx <= 10 && Nz <= 10  % save 0D plots

        name = [outpath '/',RunID,'_txc_',num2str(floor(step/nop))]; % figure 1
        print(fh1,name,'-dpng','-r300','-image');
        name = [outpath '/',RunID,'_fre_',num2str(floor(step/nop))]; % figure 2
        print(fh2,name,'-dpng','-r300','-image');

    elseif Nx <= 10   % save 1D plots

        name = [outpath '/',RunID,'_txc_',num2str(floor(step/nop))];   % figure 1
        print(fh1,name,'-dpng','-r300','-image');
        name = [outpath '/',RunID,'_fmr_',num2str(floor(step/nop))];   % figure 2
        print(fh2,name,'-dpng','-r300','-image');
        name = [outpath '/',RunID,'_vp_',num2str(floor(step/nop))];    % figure 3
        print(fh3,name,'-dpng','-r300','-image');
        name = [outpath '/',RunID,'_seg_',num2str(floor(step/nop))];   % figure 4
        print(fh4,name,'-dpng','-r300','-image');
        name = [outpath '/',RunID,'_consv_',num2str(floor(step/nop))]; % figure 5
        print(fh5,name,'-dpng','-r300','-image');

    else  % save 2D plots

        name = [outpath '/',RunID,'_vp_',num2str(floor(step/nop))]; % figure 1
        print(fh1,name,'-dpng','-r300','-image');
        name = [outpath '/',RunID,'_phsvol_',num2str(floor(step/nop))]; % figure 2
        print(fh2,name,'-dpng','-r300','-image');
        name = [outpath '/',RunID,'_phswt_',num2str(floor(step/nop))]; % figure 3
        print(fh3,name,'-dpng','-r300','-image');
        name = [outpath '/',RunID,'_vseg',num2str(floor(step/nop))]; % figure 4
        print(fh4,name,'-dpng','-r300','-image');
        name = [outpath '/',RunID,'_xc',num2str(floor(step/nop))]; % figure 5
        print(fh5,name,'-dpng','-r300','-image');

    end

    name = [outpath '/',RunID,'_phsdg',num2str(floor(step/nop))]; % figure 9 phase diagram
    print(fh9,name,'-dpng','-r300','-image');
    name = [outpath '/',RunID,'_consv',num2str(floor(step/nop))]; % figure 10 conserved quantities
    print(fh10,name,'-dpng','-r300','-image');
    name = [outpath '/',RunID,'_SDiagnostics',num2str(floor(step/nop))]; % figure 11 S and T diagnostics
    print(fh11,name,'-dpng','-r300','-image');

    name = [outpath '/',RunID,'_',num2str(step/nop)];
    save(name,'U','W','P','Pt','xFe','xSi','cFe','cSi','csFe','clFe','csSi','clSi',...
        'fsFe','flFe','fsSi','flSi','S','XFe','XSi','CFe','CSi','FsFe','FlFe','FsSi','FlSi',...
        'phisFe','philFe','phisSi','philSi','rho','Eta','segsFe','seglFe','segsSi','seglSi',...
        'Ksgr_x','Ksgr_f','Ksgr_m','xP','zP','xU','zU','xW','zW',...
        'So','XFeo','XSio','CFeo','CSio','FsFeo','FsSio',...
        'rhoo','T','yr','nxP','nzP','time','step','EntProd','dSdt','diss_T','advn_S');
    name = [outpath '/',RunID,'_cont'];
    save(name,'U','W','P','Pt','xFe','xSi','cFe','cSi','csFe','clFe','csSi','clSi',...
        'fsFe','flFe','fsSi','flSi','S','XFe','XSi','CFe','CSi','FsFe','FlFe','FsSi','FlSi',...
        'phisFe','philFe','phisSi','philSi','rho','Eta','segsFe','seglFe','segsSi','seglSi',...
        'Ksgr_x','Ksgr_f','Ksgr_m','xP','zP','xU','zU','xW','zW',...
        'So','XFeo','XSio','CFeo','CSio','FsFeo','FsSio',...
        'rhoo','T','yr','nxP','nzP','time','step','EntProd','dSdt','diss_T','advn_S');
    name = [outpath,'/',RunID,'_hist'];
    save(name,'HST');

    if step == 0
        logfile = [outpath '/',RunID,'.log'];
        if exist(logfile,'file'); delete(logfile); end
        diary(logfile)
    end
end

