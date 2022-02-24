% planetesimal: generate output

% print output header
fprintf(1,'\n*****  prepare output frame %d  *****\n',RUN.frame);


%% plot output figures
if RUN.plot
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
    imagesc(NUM.xU,NUM.zU,SOL.U); hold on;
    quiver(NUM.xP(xq),NUM.zP(zq),SOL.UP(zq,xq),SOL.WP(zq,xq),'k')
    colormap(subplot(2,3,1),cm2)
    axis ij equal tight;
    colorbar
    title('v^*_x [ms^-^1]')
    
    subplot(2,3,2)
    imagesc(NUM.xW,NUM.zW,-SOL.W); hold on;
    quiver(NUM.xP(xq),NUM.zP(zq),SOL.UP(zq,xq),SOL.WP(zq,xq),'k')
    colormap(subplot(2,3,2),cm2)
    axis ij equal tight;
    colorbar
    title('v_z^*-velocity [ms^-^1]')

    subplot(2,3,3)
    imagesc(NUM.xP,NUM.zP,SOL.P); hold on;
    colormap(subplot(2,3,3),cm2)
    axis ij equal tight;
    colorbar
    title('dynamic pressure [Pa]')
    
    Tplot = SOL.T;
    subplot(2,3,4)
    imagesc(NUM.xP,NUM.zP,SOL.T);
    colormap(subplot(2,3,4),flipud(cm1))
    axis ij equal tight;
    colorbar
    title('Temperature [C]')

    subplot(2,3,5);
    imagesc(NUM.xP,NUM.zP,MAT.rhot);
    colormap(subplot(2,3,5),cm1)
    axis ij equal tight;
    colorbar
    title('bulk density [kgm^-^3]')
    
    subplot(2,3,6);
    imagesc(NUM.xP,NUM.zP,MAT.Eta);
    colormap(subplot(2,3,6),cm1)
    axis ij equal tight;
    colorbar
    title('[\eta]')
    
%     subplot(2,3,6);
%     imagesc(NUM.xP,NUM.zP,log10(MAT.Eta));
%     colormap(subplot(2,3,6),cm1)
%     axis ij equal tight;
%     colorbar
%     title('Viscosity [log_1_0 Pas]')

% fh2 = figure(2); clf;
%     
% figure(2); imagesc(SOL.phi(end-9:end,:)); colormap(cm1); colorbar

%% plot phase vol fractions
fh2 = figure(2); clf;
figure(2) % plot phase fractions
subplot(2,2,1);
imagesc(NUM.xP,NUM.zP,SOL.phiSil);
colorbar
axis ij image
title('\phi_{Si}^l')
subplot(2,2,2);
imagesc(NUM.xP,NUM.zP,SOL.phiSis);
colorbar
axis ij image
title('\phi_{Si}^s')
subplot(2,2,3);
imagesc(NUM.xP,NUM.zP,SOL.phiFel);
colorbar
axis ij image
title('\phi_{Fe}^l')
subplot(2,2,4);
imagesc(NUM.xP,NUM.zP,SOL.phiFes);
colorbar
axis ij image
title('\phi_{Fe}^s')

%% plot phase 
fh3 = figure(3); clf
figure(3) % plot phase fractions
subplot(2,2,1);
imagesc(NUM.xP,NUM.zP,CHM.fSil);
colorbar
axis ij image
title('f_{Si}^l [wt]')
subplot(2,2,2);
imagesc(NUM.xP,NUM.zP,CHM.fSis);
colorbar
axis ij image
title('f_{Si}^s [wt]')
subplot(2,2,3);
imagesc(NUM.xP,NUM.zP,CHM.fFel);
colorbar
axis ij image
title('f_{Fe}^l [wt]')
subplot(2,2,4);
imagesc(NUM.xP,NUM.zP,CHM.fFes);
colorbar
axis ij image
title('f_{Fe}^s [wt]')

%% plot phase segregation
fh4 = figure(4); clf
figure(4) % plot phase segregation
subplot(2,2,1);
imagesc(NUM.xW,NUM.zW,SOL.phiFes+SOL.phiFel+SOL.phiSis+SOL.phiSil);
colorbar
axis ij image
title('total phi')

subplot(2,2,2);
imagesc(NUM.xP,NUM.zP,-segSis);
colorbar
axis ij image
title('v_{\Delta, Si}^s')

subplot(2,2,3);
imagesc(NUM.xP,NUM.zP,-segFel);
colorbar
axis ij image
title('v_{\Delta, Fe}^l')

subplot(2,2,4);
imagesc(NUM.xP,NUM.zP,-segFes);
colorbar
axis ij image
title('v_{\Delta, Fe}^s')

%% plot system fractions and chemical solutions
fh5 = figure(5); clf
figure(5)
subplot(2,2,1)
imagesc(NUM.xP,NUM.zP,CHM.xFe)
colorbar
axis ij image
title('x_{Fe} [wt]')

subplot(2,2,2)
imagesc(NUM.xP,NUM.zP,CHM.xSi)
colorbar
axis ij image
title('x_{Si} [wt]')

subplot(2,2,3)
imagesc(NUM.xP,NUM.zP,CHM.cFe)
colorbar
axis ij image
title('c_{Fe} [wt]')

subplot(2,2,4)
imagesc(NUM.xP,NUM.zP,CHM.cSi)
colorbar
axis ij image
title('c_{Si} [wt]')

%% plot phase densities 
fh6 = figure(6); clf
figure(6)
subplot(2,2,1)
imagesc(NUM.xP,NUM.zP,MAT.rhoSil)
colorbar
axis ij image
title('\rho_{Si}^l')

subplot(2,2,2)
imagesc(NUM.xP,NUM.zP,MAT.rhoSis)
colorbar
axis ij image
title('\rho_{Si}^s')

subplot(2,2,3)
imagesc(NUM.xP,NUM.zP,MAT.rhoFel)
colorbar
axis ij image
title('\rho_{Fe}^l')

subplot(2,2,4)
imagesc(NUM.xP,NUM.zP,MAT.rhoFes)
colorbar
axis ij image
title('\rho_{Fe}^s')

%% plot chemical phase compositions
fh7 = figure(7); clf
figure(7)
subplot(2,2,1)
imagesc(NUM.xP,NUM.zP,CHM.clSi)
colorbar
axis ij image
title('c_{Si}^l')

subplot(2,2,2)
imagesc(NUM.xP,NUM.zP,CHM.csSi)
colorbar
axis ij image
title('c_{Si}^s')

subplot(2,2,3)
imagesc(NUM.xP,NUM.zP,CHM.clFe)
colorbar
axis ij image
title('c_{Fe}^l')

subplot(2,2,4)
imagesc(NUM.xP,NUM.zP,CHM.csFe)
colorbar
axis ij image
title('c_{Fe}^s')
%% plot Partial densities and 
fh8 = figure(8); clf;
figure(8)
subplot(2,2,1)
imagesc(NUM.xP,NUM.zP,CHM.CFe)
colorbar
axis ij image
title('C_{Fe} ')

subplot(2,2,2)
imagesc(NUM.xP,NUM.zP,CHM.CSi)
colorbar
axis ij image
title('C_{Si} ')

subplot(2,2,3)
imagesc(NUM.xP,NUM.zP,CHM.XFe)
colorbar
axis ij image
title('X_{Fe}')

subplot(2,2,4)
imagesc(NUM.xP,NUM.zP,CHM.XSi)
colorbar
axis ij image
title('X_{Si}')


%% plot phase diagrams
fh9 = figure(9); clf;
TT = linspace(CHM.TSi1,CHM.TSi2,1e3);
cc = [linspace(CHM.cphsSi2,(CHM.perCsSi+CHM.perClSi)/2,round((CHM.perTSi-CHM.TSi1)./(CHM.TSi2-CHM.TSi1)*1e3)),linspace((CHM.perCsSi+CHM.perClSi)/2,CHM.cphsSi1,round((CHM.perTSi-CHM.TSi2)./(CHM.TSi1-CHM.TSi2)*1e3))];
[~,CCxSi,CClSi]     = equilibrium(TT,cc,0.*TT,CHM.TSi1,CHM.TSi2,CHM.cphsSi1,CHM.cphsSi2,...
                                  CHM.perTSi,CHM.perCsSi,CHM.perClSi,CHM.clap,CHM.PhDgSi,TINY);

TT2 = linspace(CHM.TFe1,CHM.TFe2,1000);
cc2 = linspace(CHM.cphsFe2,CHM.cphsFe1,length(TT2));
% cc2 = [linspace(cphsFe2,(CHM.perCsFe+CHM.perClFe)/2,round((CHM.perTFe-CHM.TFe1)./(CHM.TFe2-CHM.TFe1)*1e3)),linspace((CHM.perCsFe+CHM.perClFe)/2,cphsFe1,round((CHM.perTFe-CHM.TFe2)./(CHM.TFe1-CHM.TFe2)*1e3))];
[~,CCxFe,CClFe]     = equilibrium(TT2,cc2,0.*TT2,CHM.TFe1,CHM.TFe2,CHM.cphsFe1,CHM.cphsFe2,...
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
fh10 = figure(10);
% mass conservation
% if NUM.step > 0
    subplot(5,1,1)
    plot(NUM.time./3600,CON.EM(NUM.step+1),'ko','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
    ylabel('consv. $Mass$',TX{:},FS{:});
    xlabel('Time [hr]',TX{:},FS{:});
    subplot(5,1,2)
    plot(NUM.time./3600,CON.EH(NUM.step+1),'ko','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
    ylabel('consv. $H$',TX{:},FS{:}); set(gca,'XTickLabel',[]);
    subplot(5,1,3)
    plot(NUM.time./3600,CON.EX(NUM.step+1),'ko','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
    ylabel('consv. $X_{fe}$',TX{:},FS{:}); set(gca,'XTickLabel',[]);
    subplot(5,1,4)
    plot(NUM.time./3600,CON.ECFe(NUM.step+1),'ko','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
    ylabel('consv. $C_{Fe}$',TX{:},FS{:}); set(gca,'XTickLabel',[]);
    subplot(5,1,5)
    plot(NUM.time./3600,CON.ECSi(NUM.step+1),'ko','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
    ylabel('consv. $C_{Si}$',TX{:},FS{:}); set(gca,'XTickLabel',[]);
% end
    
end

%% save output
if RUN.save
    % print figure
%     name = ['../out/',RUN.ID,'/',RUN.ID,'_fig',num2str(RUN.frame)];
%     print(fh1,name,'-dpng','-r300');
   

    % save output data
    name = [outpath '/',RUN.ID,'_vp_',num2str(floor(NUM.step/RUN.nop))]; % figure 1
    print(fh1,name,'-dpng','-r300','-opengl');
    name = [outpath '/',RUN.ID,'_phsvol_',num2str(floor(NUM.step/RUN.nop))]; % figure 2
    print(fh2,name,'-dpng','-r300','-opengl');
    name = [outpath '/',RUN.ID,'_phswt_',num2str(floor(NUM.step/RUN.nop))]; % figure 3
    print(fh3,name,'-dpng','-r300','-opengl');
    name = [outpath '/',RUN.ID,'_vseg',num2str(floor(NUM.step/RUN.nop))]; % figure 4
    print(fh4,name,'-dpng','-r300','-opengl');
    name = [outpath '/',RUN.ID,'_eq',num2str(floor(NUM.step/RUN.nop))]; % figure 5
    print(fh5,name,'-dpng','-r300','-opengl');
    name = [outpath '/',RUN.ID,'_rho',num2str(floor(NUM.step/RUN.nop))]; % figure 5
    print(fh6,name,'-dpng','-r300','-opengl');
    name = [outpath '/',RUN.ID,'_chm',num2str(floor(NUM.step/RUN.nop))]; % figure 5
    print(fh7,name,'-dpng','-r300','-opengl');
    name = [outpath '/',RUN.ID,'_XC',num2str(floor(NUM.step/RUN.nop))]; % figure 5
    print(fh8,name,'-dpng','-r300','-opengl');
    name = [outpath '/',RUN.ID,'_Phsdg',num2str(floor(NUM.step/RUN.nop))]; % figure 9 phase diagram
    print(fh9,name,'-dpng','-r300','-opengl');
    name = [outpath '/',RUN.ID,'_conservation',num2str(floor(NUM.step/RUN.nop))]; % figure 10 conserved quantities
    print(fh10,name,'-dpng','-r300','-opengl');



%     name = ['../out/',RUN.ID,'/',RUN.ID,'_cont'];
%     save([name,'.mat']);
%     name = ['../out/',RUN.ID,'/',RUN.ID,'_',num2str(RUN.frame)];
%     save([name,'.mat']);

    name = [outpath '/',RUN.ID,'_',num2str(NUM.step/RUN.nop)];
    save(name,'NUM','PHY','SOL','MAT','CHM');
    name = [outpath '/',RUN.ID,'_cont'];
    save(name,'NUM','PHY','SOL','MAT','CHM');

%     % clean workspace
%     clear aa A AA advn_T diff_T dtadvn dtdiff EtaC1 EtaC2 EtaP1 EtaP2 ii II
%     clear indP indU indW IR jj JJ jj1 jj2 jj3 jj4 kappa pert Pscale RhoRef
%     clear rr R RR S To dTdto toc_assmb toc_solve toc_update xq zq V fh1
    
    if NUM.step == 0
        logfile = [outpath '/',RUN.ID,'.log'];
        if exist(logfile,'file'); delete(logfile); end
        diary(logfile)
    end
end

RUN.frame = RUN.frame + 1;  % increment frame count
