% planetesimal: generate output

% print output header
fprintf(1,'*****  prepare output frame %d  *****\n',RUN.frame);


%% plot output figures
if RUN.plot
    
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
imagesc(NUM.xP,NUM.zP,phiSil);
colorbar
axis ij image
title('\phi_{Si}^l')
subplot(2,2,2);
imagesc(NUM.xP,NUM.zP,phiSis);
colorbar
axis ij image
title('\phi_{Si}^s')
subplot(2,2,3);
imagesc(NUM.xP,NUM.zP,phiFel);
colorbar
axis ij image
title('\phi_{Fe}^l')
subplot(2,2,4);
imagesc(NUM.xP,NUM.zP,phiFes);
colorbar
axis ij image
title('\phi_{Fe}^s')

%% plot phase 
fh3 = figure(3); clf
figure(3) % plot phase fractions
subplot(2,2,1);
imagesc(NUM.xP,NUM.zP,fSil);
colorbar
axis ij image
title('f_{Si}^l [wt]')
subplot(2,2,2);
imagesc(NUM.xP,NUM.zP,fSis);
colorbar
axis ij image
title('f_{Si}^s [wt]')
subplot(2,2,3);
imagesc(NUM.xP,NUM.zP,fFel);
colorbar
axis ij image
title('f_{Fe}^l [wt]')
subplot(2,2,4);
imagesc(NUM.xP,NUM.zP,fFes);
colorbar
axis ij image
title('f_{Fe}^s [wt]')

%% plot phase segregation
fh4 = figure(4); clf
figure(4) % plot phase segregation
subplot(2,2,1);
imagesc(NUM.xW,NUM.zW,phiFes+phiFel+phiSis+phiSil);
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
imagesc(NUM.xP,NUM.zP,xFe)
colorbar
axis ij image
title('x_{Fe} [wt]')

subplot(2,2,2)
imagesc(NUM.xP,NUM.zP,xSi)
colorbar
axis ij image
title('x_{Si} [wt]')

subplot(2,2,3)
imagesc(NUM.xP,NUM.zP,cFe)
colorbar
axis ij image
title('c_{Fe} [wt]')

subplot(2,2,4)
imagesc(NUM.xP,NUM.zP,cSi)
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
imagesc(NUM.xP,NUM.zP,clSi)
colorbar
axis ij image
title('c_{Si}^l')

subplot(2,2,2)
imagesc(NUM.xP,NUM.zP,csSi)
colorbar
axis ij image
title('c_{Si}^s')

subplot(2,2,3)
imagesc(NUM.xP,NUM.zP,clFe)
colorbar
axis ij image
title('c_{Fe}^l')

subplot(2,2,4)
imagesc(NUM.xP,NUM.zP,csFe)
colorbar
axis ij image
title('c_{Fe}^s')
%% plot Partial densities and 
fh8 = figure(8); clf;
figure(8)
subplot(2,2,1)
imagesc(NUM.xP,NUM.zP,CFe)
colorbar
axis ij image
title('C_{Fe} ')

subplot(2,2,2)
imagesc(NUM.xP,NUM.zP,CSi)
colorbar
axis ij image
title('C_{Si} ')

subplot(2,2,3)
imagesc(NUM.xP,NUM.zP,XFe)
colorbar
axis ij image
title('X_{Fe}')

subplot(2,2,4)
imagesc(NUM.xP,NUM.zP,XSi)
colorbar
axis ij image
title('X_{Si}')


%% plot phase diagrams
fh9 = figure(9); clf;
TT = linspace(TSi1,TSi2,1e3);
cc = [linspace(cphsSi2,(perCsSi+perClSi)/2,round((perTSi-TSi1)./(TSi2-TSi1)*1e3)),linspace((perCsSi+perClSi)/2,cphsSi1,round((perTSi-TSi2)./(TSi1-TSi2)*1e3))];
[~,CCxSi,CClSi]     = equilibrium(TT,cc,0.*TT,TSi1,TSi2,cphsSi1,cphsSi2,...
                                  perTSi,perCsSi,perClSi,clap,PhDgSi,TINY);

TT2 = linspace(TFe1,TFe2,1000);
cc2 = linspace(cphsFe2,cphsFe1,length(TT2));
% cc2 = [linspace(cphsFe2,(perCsFe+perClFe)/2,round((perTFe-TFe1)./(TFe2-TFe1)*1e3)),linspace((perCsFe+perClFe)/2,cphsFe1,round((perTFe-TFe2)./(TFe1-TFe2)*1e3))];
[~,CCxFe,CClFe]     = equilibrium(TT2,cc2,0.*TT2,TFe1,TFe2,cphsFe1,cphsFe2,...
                                  perTFe,perCsFe,perClFe,clap,PhDgFe,TINY);
subplot(1,2,1)
plot(CCxSi,TT,'k-','LineWidth',2); axis tight; hold on; box on;
plot(CClSi,TT,'k-','LineWidth',2); axis tight; hold on; box on;
plot(cSi(2:end-1),SOL.T(2:end-1), '.k', csSi(2:end-1),SOL.T(2:end-1), '.b', clSi(2:end-1),SOL.T(2:end-1), '.r','MarkerSize',15)


subplot(1,2,2)
plot(CCxFe,TT2,'k-','LineWidth',2); axis tight; hold on; box on;
plot(CClFe,TT2,'k-','LineWidth',2); axis tight; hold on; box on;
plot(cFe(2:end-1),SOL.T(2:end-1), '.k', csFe(2:end-1),SOL.T(2:end-1), '.b', clFe(2:end-1),SOL.T(2:end-1), '.r','MarkerSize',15)

drawnow;

    
    
end

%% save output
if RUN.save
    % print figure
    name = ['../out/',RUN.ID,'/',RUN.ID,'_fig',num2str(RUN.frame)];
    print(fh1,name,'-dpng','-r300');
    
    % clean workspace
    clear aa A AA advn_T diff_T dtadvn dtdiff EtaC1 EtaC2 EtaP1 EtaP2 ii II
    clear indP indU indW IR jj JJ jj1 jj2 jj3 jj4 kappa pert Pscale RhoRef
    clear rr R RR S To dTdto toc_assmb toc_solve toc_update xq zq V fh1

    % save output data
    name = ['../out/',RUN.ID,'/',RUN.ID,'_cont'];
    save([name,'.mat']);
    name = ['../out/',RUN.ID,'/',RUN.ID,'_',num2str(RUN.frame)];
    save([name,'.mat']);
    
    if NUM.step == 0
        logfile = ['../out/',RUN.ID,'/',RUN.ID,'.log'];
        if exist(logfile,'file'); delete(logfile); end
        diary(logfile)
    end
end

RUN.frame = RUN.frame + 1;  % increment frame count
