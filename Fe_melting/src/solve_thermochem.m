% planetesimal: solve thermo-chemical equations


%% update enthalpy and Temperature
advn_H  = advection(MAT.rho.*CHM.xSi.*CHM.fSis.*SOL.T.* PHY.CpSis               ,UsSi,WsSi,NUM.h,NUM.h,NUM.ADVN,'flx')...
        + advection(MAT.rho.*CHM.xFe.*CHM.fFes.*SOL.T.* PHY.CpFes               ,UsFe,WsFe,NUM.h,NUM.h,NUM.ADVN,'flx')...
        + advection(MAT.rho.*CHM.xSi.*CHM.fSil.*SOL.T.*(PHY.CpSil + CHM.dEntrSi),UlSi,WlSi,NUM.h,NUM.h,NUM.ADVN,'flx')...
        + advection(MAT.rho.*CHM.xFe.*CHM.fFel.*SOL.T.*(PHY.CpFel + CHM.dEntrFe),UlFe,WlFe,NUM.h,NUM.h,NUM.ADVN,'flx');

qTz    = - (MAT.kT(1:end-1,:)+MAT.kT(2:end,:))./2 .* ddz(SOL.T,NUM.h);     % heat diffusion z-flux
qTx    = - (MAT.kT(:,1:end-1)+MAT.kT(:,2:end))./2 .* ddx(SOL.T,NUM.h);     % heat diffusion x-flux
diff_T(2:end-1,2:end-1) = (- ddz(qTz(:,2:end-1),NUM.h) ...                 % heat diffusion
                           - ddx(qTx(2:end-1,:),NUM.h));

dHdt   = - advn_H + diff_T + MAT.Hr;

if NUM.step>0
    % update solution
    SOL.H = Ho + (NUM.theta.*dHdt + (1-NUM.theta).*dHdto) .* NUM.dt;
    
    %apply boundaries
    switch SOL.BCTTop
        case 'isothermal'
            SOL.H(1,:)  =  SOL.T0.*(MAT.rhoDs(1,:) + MAT.rhoCp(1,:));
        case 'insulating'
            SOL.H(1,:)  =  SOL.H(2,:);
    end
    switch SOL.BCTBot
        case 'isothermal'
            SOL.H(end,:)  =  SOL.T0.*(MAT.rhoDs(end,:) + MAT.rhoCp(end,:));
        case 'insulating'
            SOL.H(end,:)  =  SOL.H(end-1,:);
    end    
    switch SOL.BCTSides
        case 'isothermal'
            SOL.H(:,[1 end])  =  SOL.T0.*(MAT.rhoDs(:,[1 end]) + MAT.rhoCp(:,[1 end]));
        case 'insulating'
            SOL.H(:,[1 end])  =  SOL.H(:,[2 end-1]);
    end 
    %visualise changes
    %     figure(30)
    %     subplot(2,2,1)
    %     imagesc(Ho); colorbar; axis ij equal tight; title('Ho')
    %     subplot(2,2,2)
    %     imagesc(SOL.H); colorbar; axis ij equal tight; title('H')
    %     subplot(2,2,3)
    %     imagesc(SOL.H-Ho); colorbar; axis ij equal tight; title('dH')
    %     subplot(2,2,4)
    %     imagesc(dHdt); colorbar; axis ij equal tight; title('dH/dt')
    %     drawnow
end


%% update composition

% update system, only one system needs to be solved as SUM_i(X_i) = 1
% equation 7
adv_XFe     = advection(MAT.rho.*CHM.xFe.*CHM.fFes,UsFe,WsFe,NUM.h,NUM.h,NUM.ADVN,'flx')...
            + advection(MAT.rho.*CHM.xFe.*CHM.fFel,UlFe,WlFe,NUM.h,NUM.h,NUM.ADVN,'flx');
dXdt        = - adv_XFe;

if NUM.step>0
    % update solution
    CHM.XFe = XFeo + (NUM.theta.*dXdt + (1-NUM.theta).*dXdto) .* NUM.dt;
    
    % apply boundaries
    CHM.XFe([1 end],:) = CHM.XFe([2 end-1],:);  CHM.XFe(:,[1 end]) = CHM.XFe(:,[2 end-1]);
    
%     figure(40)
%     subplot(2,2,1); imagesc(XFeo); colorbar; title('Xo')
%     subplot(2,2,2); imagesc(CHM.XFe); colorbar; title('X')
%     subplot(2,2,3); imagesc(CHM.XFe - XFeo); colorbar; title('dX')
%     subplot(2,2,4); imagesc(dXdt); colorbar; title('dXdt')
end

% update fertile chemical components
% equations 8a, 8b
advn_CSi  = advection(MAT.rho.*CHM.xSi.*CHM.fSis.*CHM.csSi,UsSi,WsSi,NUM.h,NUM.h,NUM.ADVN,'flx')...
          + advection(MAT.rho.*CHM.xSi.*CHM.fSil.*CHM.clSi,UlSi,WlSi,NUM.h,NUM.h,NUM.ADVN,'flx');

qCSiz     = - PHY.kC.* (MAT.rho (1:end-1,:)+MAT.rho (2:end,:))/2 ...
                    .* (CHM.xSi (1:end-1,:)+CHM.xSi (2:end,:))/2 ...
                    .* (CHM.fSil(1:end-1,:)+CHM.fSil(2:end,:))/2 ...
                    .* ddz(CHM.cSi,NUM.h);                                 % chemical diffusion z-flux
qCSix     = - PHY.kC.* (MAT.rho (:,1:end-1)+MAT.rho (:,2:end))/2 ...
                    .* (CHM.xSi (:,1:end-1)+CHM.xSi (:,2:end))/2 ...
                    .* (CHM.fSil(:,1:end-1)+CHM.fSil(:,2:end))/2 ...
                    .* ddx(CHM.cSi,NUM.h); 
diff_CSi(2:end-1,2:end-1) = (- ddz(qCSiz(:,2:end-1),NUM.h) ...             % chemical diffusion
                             - ddx(qCSix(2:end-1,:),NUM.h));

dCSidt    = - advn_CSi + diff_CSi;

advn_CFe  = advection(MAT.rho.*CHM.xFe.*CHM.fFes.*CHM.csFe,UsFe,WsFe,NUM.h,NUM.h,NUM.ADVN,'flx')...
          + advection(MAT.rho.*CHM.xFe.*CHM.fFel.*CHM.clFe,UlFe,WlFe,NUM.h,NUM.h,NUM.ADVN,'flx');

qCFez     = - PHY.kC.* (MAT.rho (1:end-1,:)+MAT.rho (2:end,:))/2 ...
                    .* (CHM.xFe (1:end-1,:)+CHM.xFe (2:end,:))/2 ...
                    .* (CHM.fFel(1:end-1,:)+CHM.fFel(2:end,:))/2 ...
                    .* ddz(CHM.cFe,NUM.h);                                 % chemical diffusion z-flux
qCFex     = - PHY.kC.* (MAT.rho (:,1:end-1)+MAT.rho (:,2:end))/2 ...
                    .* (CHM.xFe (:,1:end-1)+CHM.xFe (:,2:end))/2 ...
                    .* (CHM.fFel(:,1:end-1)+CHM.fFel(:,2:end))/2 ...
                    .* ddx(CHM.cFe,NUM.h);                                 % chemical diffusion z-flux
diff_CFe(2:end-1,2:end-1) = (- ddz(qCFez(:,2:end-1),NUM.h) ...             % chemical diffusion
                             - ddx(qCFex(2:end-1,:),NUM.h));

dCFedt      = - advn_CFe + diff_CFe;

if NUM.step>0
    % update solution
    CHM.CSi = CSio + (NUM.theta.*dCSidt + (1-NUM.theta).*dCSidto) .* NUM.dt;
    CHM.CFe = CFeo + (NUM.theta.*dCFedt + (1-NUM.theta).*dCFedto) .* NUM.dt;
    
    % apply boundaries
    CHM.CSi([1 end],:) = CHM.CSi([2 end-1],:);  CHM.CSi(:,[1 end]) = CHM.CSi(:,[2 end-1]);
    CHM.CFe([1 end],:) = CHM.CFe([2 end-1],:);  CHM.CFe(:,[1 end]) = CHM.CFe(:,[2 end-1]);

%     figure(41)
%     subplot(2,2,1); imagesc(CFeo); colorbar; title('CFeo')
%     subplot(2,2,2); imagesc(CHM.CFe); colorbar; title('CFe')
%     subplot(2,2,3); imagesc(CHM.CFe - CFeo); colorbar; title('dCFe')
%     subplot(2,2,4); imagesc(dCFedt); colorbar; title('dCFedt')
% 
%     figure(42)
%     subplot(2,2,1); imagesc(CSio); colorbar; title('CSio')
%     subplot(2,2,2); imagesc(CHM.CSi); colorbar; title('CSi')
%     subplot(2,2,3); imagesc(CHM.CSi - CSio); colorbar; title('dCSi')
%     subplot(2,2,4); imagesc(dCSidt); colorbar; title('dCSidt')
%     drawnow
end

if NUM.step>0
    CHM.XSi = MAT.rho - CHM.XFe;
    CHM.xFe = CHM.XFe./MAT.rho;    
    CHM.xSi = 1-CHM.xFe;
    CHM.cSi = CHM.CSi./(CHM.xSi+TINY)./MAT.rho;
    CHM.cFe = CHM.CFe./(CHM.xFe+TINY)./MAT.rho;
    SOL.T   = SOL.H./(MAT.rhoCp + MAT.rhoDs);
else
%     SOL.H   = SOL.T.*(MAT.rhoDs + MAT.rhoCp);                              % mixture enthalpy density
%     CHM.XFe = MAT.rho.*CHM.xFe; CHM.XSi = MAT.rho.*CHM.xSi;                % mixture Fe/Si system densities
%     CHM.CFe = MAT.rho.*CHM.cFe.*CHM.xFe;                                   % mixture Fe component density
%     CHM.CSi = MAT.rho.*CHM.cSi.*CHM.xSi;                                   % mixture Si component density
end


%% update local phase equilibrium
if NUM.step>0
    % [CHM.fFes,CHM.csFe,CHM.clFe]     = equilibrium((SOL.T+To)./2,(CHM.cFe+cFeo)./2,(SOL.Pt+Pto)./2,CHM.TFe1,CHM.TFe2,CHM.cphsFe1,CHM.cphsFe2,...
    %     CHM.perTFe,CHM.perCsFe,CHM.perClFe,CHM.clap,CHM.PhDgFe,TINY);
    %
    % [CHM.fSis,CHM.csSi,CHM.clSi]     = equilibrium((SOL.T+To)./2,(CHM.cSi+cSio)./2,(SOL.Pt+Pto)./2,CHM.TSi1,CHM.TSi2,CHM.cphsSi1,CHM.cphsSi2,...
    %     CHM.perTSi,CHM.perCsSi,CHM.perClSi,CHM.clap,CHM.PhDgSi,TINY);
    % CHM.fFel = 1-CHM.fFes; CHM.fSil = 1-CHM.fSis;
    
    [fFesq,csFeq,clFeq] = equilibrium(SOL.T,CHM.cFe,SOL.Pt,CHM.TFe1,CHM.TFe2,CHM.cphsFe1,CHM.cphsFe2,...
                                      CHM.perTFe,CHM.perCsFe,CHM.perClFe,CHM.clap,CHM.PhDgFe,TINY);
    
    [fSisq,csSiq,clSiq] = equilibrium(SOL.T,CHM.cSi,SOL.Pt,CHM.TSi1,CHM.TSi2,CHM.cphsSi1,CHM.cphsSi2,...
                                      CHM.perTSi,CHM.perCsSi,CHM.perClSi,CHM.clap,CHM.PhDgSi,TINY);
    
    % lag equilibrium phase fractions
    CHM.fFes = alpha.*CHM.fFes + (1-alpha).*fFesq;
    CHM.fSis = alpha.*CHM.fSis + (1-alpha).*fSisq;
    
    CHM.fFel = 1-CHM.fFes; CHM.fSil = 1-CHM.fSis;
    
    % update phase compositions
    KcFe = csFeq./clFeq;
    CHM.clFe = CHM.cFe./(CHM.fFel + CHM.fFes.*KcFe);
    CHM.csFe = CHM.cFe./(CHM.fFel./KcFe + CHM.fFes);
    
    KcSi = csSiq./clSiq;
    CHM.clSi = CHM.cSi./(CHM.fSil + CHM.fSis.*KcSi);
    CHM.csSi = CHM.cSi./(CHM.fSil./KcSi + CHM.fSis);
    
end

% %check outputs
% NUM.dt
%
% figure(101)
% subplot(2,2,1)
% imagesc(SOL.H-Ho); colorbar
% title('dH')
%
% subplot(2,2,2)
% imagesc(CHM.XFe - CHM.XFeo); colorbar
% title('dX')
%
% subplot(2,2,3)
% imagesc(CHM.CSi - CSio); colorbar
% title('dCSi')
%
% subplot(2,2,4)
% imagesc(CHM.CFe - CFeo); colorbar
% title('dCFe')
%
% drawnow

%%

% %% update phase fractions
% WFe                 = SOL.segW + SOL.W;
% [advn]              = advection(SOL.phi,SOL.U,WFe,NUM.h,NUM.h,NUM.ADVN,'Flx');
% if NUM.step>0
% SOL.phi             = phio - advn * NUM.dt;
% end
% SOL.phi(SOL.phi<0)  = 0;
% %% Solve energy equation explicitly
%
% % get advection rate
% advn_T              = advection(SOL.T,SOL.U,SOL.W,NUM.h,NUM.h,NUM.ADVN,'adv');
% % get diffusion rate
% diff_T              = (diff(SOL.T(:,2:end-1),2,1)./NUM.h^2  ...
%                     +  diff(SOL.T(2:end-1,:),2,2)./NUM.h^2) ...
%                     .* MAT.kT(2:end-1,2:end-1);
%
% % heat capacity density
% RhoCp               = MAT.Rhot.*MAT.Cp;
%
% % get total rate of change
% SOL.dTdt            = - advn_T(2:end-1,2:end-1) ...
%                       + (diff_T + MAT.Hr(2:end-1,2:end-1) ...
%                       + SOL.Hs(2:end-1,2:end-1) ...
%                       + SOL.Ha(2:end-1,2:end-1))./RhoCp(2:end-1,2:end-1);
% SOL.dTdt(isinf(SOL.dTdt)) = 0;
%
%
% % update temperature solution
% if NUM.step>0
% SOL.T(2:end-1,2:end-1) = To(2:end-1,2:end-1) + (  NUM.theta .*SOL.dTdt   ...
%     + (1-NUM.theta).*    dTdto) .* NUM.dt;
% end
%
% % apply top boundary conditions
% switch SOL.BCTTop
%     case 'isothermal'; SOL.T(1,:) =     To(1,:);
%     case 'insulating'; SOL.T(1,:) = SOL.T (2,:);
% end
%
% % apply bottom boundary conditions
% switch SOL.BCTBot
%     case 'isothermal'; SOL.T(end,:) =     To(end  ,:);
%     case 'insulating'; SOL.T(end,:) = SOL.T (end-1,:);
% end
%
% % apply side boundary conditions
% switch SOL.BCTSides
%     case 'isothermal'; SOL.T(:,[1 end]) =     To(:,[1 end  ]);
%     case 'insulating'; SOL.T(:,[1 end]) = SOL.T (:,[2 end-1]);
% end