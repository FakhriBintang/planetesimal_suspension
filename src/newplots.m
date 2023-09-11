%% newplots
% a collection of post-processing plots separate from the default plotsfor better visualisation in
% manuscripts
% add path to source directory
addpath('../src')
addpath('../src/cbrewer/')

% use color brewer to create colormaps
cm1 =        cbrewer('seq','YlOrRd',30) ; % sequential colour map
cm2 = flipud(cbrewer('div','RdBu'  ,30)); % divergent colour map

TX = {'Interpreter','Latex'}; FS = {'FontSize',16};
TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',14};
UN = {'Units','Centimeters'};
CL = {'Color',[0.0 0.0 0.0],[0.80 0.15 0.10],[0.10 0.15 0.65],[0.45 0.60 0.95]};
LW = {'LineWidth',2};

%% row averaging plots
% to show more difference between

Fe_avg = mean(rho(2:end-1,2:end-1),2);
subplot(2,2,1); imagesc(xP(2:end-1),zP(2:end-1),rho(2:end-1,2:end-1)-Fe_avg)
colorbar
axis ij equal tight
title('lateral rho [wt] rho- row averaged rho')
subplot(2,2,3); imagesc(xP(2:end-1),zP(2:end-1),rho(2:end-1,2:end-1))
colorbar
axis ij equal tight
title('true rho [wt]')
subplot(2,2,[2,4]);plot(mean(rho(2:end-1,2:end-1),2),zP(2:end-1),'LineWidth',2); axis ij tight; box on;
title('row averaged rho [wt]')

% %% temporal evolution per row
% % paste one row at a time
% subplot(3,3,1)
% imagesc(xP(2:end-1),zP(2:end-1),T(2:end-1,2:end-1));
% colormap(subplot(3,3,1),flipud(cm1))
% axis ij equal tight;
% colorbar
% title('Temperature [C]',TX{:},FS{:}); set(gca,TL{:},TS{:});
% 
% subplot(3,3,2);
% imagesc(xP(2:end-1),zP(2:end-1),rho(2:end-1,2:end-1));
% colormap(subplot(3,3,2),cm1)
% axis ij equal tight;
% colorbar
% title('bulk density $\bar{\rho}$ $[kgm^{-3}$]',TX{:},FS{:}); set(gca,TL{:},TS{:});
% 
% subplot(3,3,3);
% imagesc(xP(2:end-1),zP(2:end-1),log10(Eta(2:end-1,2:end-1)));
% colormap(subplot(3,3,3),cm1)
% axis ij equal tight;
% colorbar
% title('\eta $[log_{10} Pas]$',TX{:},FS{:}); set(gca,TL{:},TS{:});
% % row 2
% subplot(3,3,4)
% imagesc(xP(2:end-1),zP(2:end-1),T(2:end-1,2:end-1));
% colormap(subplot(3,3,4),flipud(cm1))
% axis ij equal tight;
% colorbar
% title('Temperature [C]',TX{:},FS{:}); set(gca,TL{:},TS{:});
% 
% subplot(3,3,5);
% imagesc(xP(2:end-1),zP(2:end-1),rho(2:end-1,2:end-1));
% colormap(subplot(3,3,5),cm1)
% axis ij equal tight;
% colorbar
% title('bulk density $\bar{\rho}$ $[kgm^{-3}$]',TX{:},FS{:}); set(gca,TL{:},TS{:});
% 
% subplot(3,3,6);
% imagesc(xP(2:end-1),zP(2:end-1),log10(Eta(2:end-1,2:end-1)));
% colormap(subplot(3,3,6),cm1)
% axis ij equal tight;
% colorbar
% title('$\eta$ $[log_{10} Pas]$',TX{:},FS{:}); set(gca,TL{:},TS{:});
% % row 3
% subplot(3,3,7)
% imagesc(xP(2:end-1),zP(2:end-1),T(2:end-1,2:end-1));
% colormap(subplot(3,3,7),flipud(cm1))
% axis ij equal tight;
% colorbar
% title('Temperature [C]',TX{:},FS{:}); set(gca,TL{:},TS{:});
% 
% subplot(3,3,8);
% imagesc(xP(2:end-1),zP(2:end-1),rho(2:end-1,2:end-1));
% colormap(subplot(3,3,8),cm1)
% axis ij equal tight;
% colorbar
% title('bulk density $\bar{\rho}$ $[kgm^{-3}$]',TX{:},FS{:}); set(gca,TL{:},TS{:});
% 
% subplot(3,3,9);
% imagesc(xP(2:end-1),zP(2:end-1),log10(Eta(2:end-1,2:end-1)));
% colormap(subplot(3,3,9),cm1)
% axis ij equal tight;
% colorbar
% title('$\eta$ $[log_{10} Pas]$',TX{:},FS{:}); set(gca,TL{:},TS{:});