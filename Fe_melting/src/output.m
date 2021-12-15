% planetesimal: generate output

% print output header
fprintf(1,'*****  prepare output frame %d  *****\n',RUN.frame);


%% plot output figures
if RUN.plot
    
    xq = round(NUM.N/10/2)+1:round(NUM.N/10):NUM.nxP;  % x-indexes for quiver plots
    zq = round(NUM.N/10/2)+1:round(NUM.N/10):NUM.nxP;  % z-indexes for quiver plots
    
    fh1 = figure(1); clf;
    
    fh2 = figure(2); clf; colormap(cm2)
    
    fh3 = figure(3); clf;
    
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
    
%     subplot(2,3,3)
%     imagesc(NUM.xP,NUM.zP,log10(SOL.Stokes)); hold on;
% %     imagesc(NUM.xP,NUM.zP,SOL.Stokes); hold on;
%     quiver(NUM.xP(xq),NUM.zP(zq),SOL.UP(zq,xq),SOL.WP(zq,xq),'k')
%     colormap(subplot(2,3,3),cm2)
%     axis ij equal tight;
%     colorbar
%     title('log particle settling velocity')

    subplot(2,3,3)
    imagesc(NUM.xP,NUM.zP,SOL.P); hold on;
%     quiver(NUM.xP(xq),NUM.zP(zq),SOL.UP(zq,xq),SOL.WP(zq,xq),'k')
    colormap(subplot(2,3,3),cm2)
    axis ij equal tight;
    colorbar
    title('dynamic pressure [Pa]')
    
    Tplot = SOL.T;
%     Tplot(NUM.PHI<=0) = NaN;
    subplot(2,3,4)
    imagesc(NUM.xP,NUM.zP,SOL.T);
    colormap(subplot(2,3,4),flipud(cm1))
    axis ij equal tight;
    colorbar
    title('Temperature [C]')
    
%     Rhoplot = MAT.Rho;
%     Rhoplot(NUM.PHI<=0) = NaN;
    subplot(2,3,5);
    imagesc(NUM.xP,NUM.zP,MAT.rhot);
    colormap(subplot(2,3,5),cm1)
    axis ij equal tight;
    colorbar
    title('bulk density [kgm^-^3]')
    
%     subplot(2,3,6);
%     imagesc(NUM.xP,NUM.zP,MAT.rhoSil.*phiSil + MAT.rhoSis.*phiSis + MAT.rhoFel.*phiFel +MAT.rhoFes.*phiFes);
%     colormap(subplot(2,3,6),cm1)
%     axis ij equal tight;
%     colorbar
%     title('bulk density phi method [kgm^-^3]')
    
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

figure(3) % plot phase segregation
subplot(2,2,1);
imagesc(NUM.xW,NUM.zW,phiFes+phiFel+phiSis+phiSil);
colorbar
axis ij image
title('total phi')

subplot(2,2,2);
imagesc(NUM.xP,NUM.zP,segSis);
colorbar
axis ij image
title('v_{\Delta, Si}^s')

subplot(2,2,3);
imagesc(NUM.xP,NUM.zP,segFel);
colorbar
axis ij image
title('v_{\Delta, Fe}^l')

subplot(2,2,4);
imagesc(NUM.xP,NUM.zP,segFes);
colorbar
axis ij image
title('v_{\Delta, Fe}^s')

figure(4)
subplot(2,2,1)
imagesc(NUM.xP,NUM.zP,xFe)
colorbar
axis ij image
title('Fe wt fraction')

subplot(2,2,2)
imagesc(NUM.xP,NUM.zP,xSi)
colorbar
axis ij image
title('Si wt fraction')

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

figure(5)
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

figure(6)
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

    drawnow;

    
    
end

% save output
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
