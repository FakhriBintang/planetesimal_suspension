% planetesimal: solve thermo-chemical equations

% print solver header
fprintf(1,'  ---  solve thermo-chemical equations \n');
tic;

% bndmode = 0;
%
% switch bndmode
%     case 0
%         bndshape = zeros(size(SOL.T));
%     case 1
%         bndshape = exp( ( -ZZ)/dw);
%     case 2
%         bndshape = exp( ( -ZZ)/dw) ...
%                  + exp(-(D-ZZ)/dw);
%     case 3
%         bndshape = exp( ( -ZZ)/dw) ...
%                  + exp(-(D-ZZ)/dw) ...
%                  + exp( ( -XX)/dw) ...
%                  + exp(-(L-XX)/dw);
% end


%% update enthalpy and Temperature
advn_H      = advection(MAT.rhot.*CHM.xSi.*CHM.fSis.*PHY.CpSis.*SOL.T,...
                        UsSi,WsSi,NUM.h,NUM.h,NUM.AdvnScheme,'flx')...
            + advection(MAT.rhot.*CHM.xFe.*CHM.fFes.*PHY.CpFes.*SOL.T,...
                        UsFe,WsFe,NUM.h,NUM.h,NUM.AdvnScheme,'flx')...
            + advection(MAT.rhot.*CHM.xSi.*CHM.fSil.*SOL.T.*(PHY.CpSil + CHM.dEntrSi),...
                        UlSi,WlSi,NUM.h,NUM.h,NUM.AdvnScheme,'flx')...
            + advection(MAT.rhot.*CHM.xFe.*CHM.fFel.*SOL.T.*(PHY.CpFel + CHM.dEntrFe),...
                        UlFe,WlFe,NUM.h,NUM.h,NUM.AdvnScheme,'flx');

qTz    = - (MAT.kT(1:end-1,:)+MAT.kT(2:end,:))./2 .* ddz(SOL.T,NUM.h);                 % heat diffusion z-flux
qTx    = - (MAT.kT(:,1:end-1)+MAT.kT(:,2:end))./2 .* ddx(SOL.T,NUM.h);                 % heat diffusion x-flux
diff_T = zeros(NUM.N+2,NUM.N+2);
diff_T(2:end-1,2:end-1) = (- ddz(qTz(:,2:end-1),NUM.h) ...                 % heat diffusion
    - ddx(qTx(2:end-1,:),NUM.h));
dHdt        = diff_T - advn_H +MAT.Hr +zeros(size(SOL.T));

if NUM.step>0
    SOL.H = Ho + (NUM.theta.*dHdt   + (1-NUM.theta).*dHdto)  .*NUM.dt;
    %apply boundaries

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
    switch Topbound
        case 'constant_T'
        rhoCpt([1 ],:)  = (MAT.rhot([1 ],:).*CHM.xFe([1 ],:).*(CHM.fFes([1 ],:).*PHY.CpFes + CHM.fFel([1 ],:).*PHY.CpFel)...
            +  MAT.rhot([1 ],:).*CHM.xSi([1 ],:).*(CHM.fSis([1 ],:).*PHY.CpSis + CHM.fSil([1 ],:).*PHY.CpSil));
        SOL.H([1 ],:)   =  SOL.T0.*(MAT.rhot([1 ],:).*CHM.xFe([1 ],:).*(CHM.fFel([1 ],:).*CHM.dEntrSi) ...
            +  MAT.rhot([1 ],:).*CHM.xSi([1 ],:).*(CHM.fSil([1 ],:).*CHM.dEntrSi) +rhoCpt([1 ],:));

        otherwise
        SOL.H([1 end],:) = SOL.H([2 end-1],:);  SOL.H(:,[1 end]) = SOL.H(:,[2 end-1]);
    end
end

%% update composition

% update system, only one system needs to be solved as SUM_i(X_i) = 1
% equation 7
adv_XFe     = advection(MAT.rhot.*CHM.xFe.*CHM.fFes,UsFe,WsFe,NUM.h,NUM.h,NUM.AdvnScheme,'flx')...
    + advection(MAT.rhot.*CHM.xFe.*CHM.fFel,UlFe,WlFe,NUM.h,NUM.h,NUM.AdvnScheme,'flx');
dXdt        = -adv_XFe;

if NUM.step>0
    CHM.XFe = XFeo + (NUM.theta.*dXdt   + (1-NUM.theta).*dXdto)  .*NUM.dt;
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
advn_CSi    = advection(MAT.rhot.*CHM.xSi.*CHM.fSis.*CHM.csSi,...
    UsSi,WsSi,NUM.h,NUM.h,NUM.AdvnScheme,'flx')...
    + advection(MAT.rhot.*CHM.xSi.*CHM.fSil.*CHM.clSi,...
    UlSi,WlSi,NUM.h,NUM.h,NUM.AdvnScheme,'flx');

qCSiz    = - (MAT.kC(1:end-1,:).*CHM.fSil(1:end-1,:).*MAT.rhot(1:end-1,:) ...
    +MAT.kC(2:end,:)  .*CHM.fSil(2:end,:)  .*MAT.rhot(2:end,:))./2 ...
    .* ddz(CHM.clSi,NUM.h);                 % heat diffusion z-flux
qCSix    = - (MAT.kC(:,1:end-1).*CHM.fSil(:,1:end-1).*MAT.rhot(:,1:end-1)...
    + MAT.kC(:,2:end)  .*CHM.fSil(:,2:end)  .*MAT.rhot(:,2:end))./2 ...
    .* ddx(CHM.clSi,NUM.h);                 % heat diffusion x-flux
diff_CSi(2:end-1,2:end-1) = (- ddz(qCSiz(:,2:end-1),NUM.h) ...                 % heat diffusion
    - ddx(qCSix(2:end-1,:),NUM.h));

dCSidt      = -advn_CSi + diff_CSi +zeros(size(CHM.CSi));

advn_CFe    = advection(MAT.rhot.*CHM.xFe.*CHM.fFes.*CHM.csFe,...
    UsFe,WsFe,NUM.h,NUM.h,NUM.AdvnScheme,'flx')...
    + advection(MAT.rhot.*CHM.xFe.*CHM.fFel.*CHM.clFe,...
    UlFe,WlFe,NUM.h,NUM.h,NUM.AdvnScheme,'flx');

qCFez    = - (MAT.kC(1:end-1,:).*CHM.fFel(1:end-1,:).*MAT.rhot(1:end-1,:) ...
    +MAT.kC(2:end,:)  .*CHM.fFel(2:end,:)  .*MAT.rhot(2:end,:))./2 ...
    .* ddz(CHM.clFe,NUM.h);                 % heat diffusion z-flux
qCFex    = - (MAT.kC(:,1:end-1).*CHM.fFel(:,1:end-1).*MAT.rhot(:,1:end-1)...
    + MAT.kC(:,2:end)  .*CHM.fFel(:,2:end)  .*MAT.rhot(:,2:end))./2 ...
    .* ddx(CHM.clFe,NUM.h);                 % heat diffusion x-flux
diff_CFe(2:end-1,2:end-1) = (- ddz(qCFez(:,2:end-1),NUM.h) ...                 % heat diffusion
    - ddx(qCFex(2:end-1,:),NUM.h));

dCFedt      = -advn_CFe + diff_CFe + zeros(size(CHM.CFe));

if NUM.step>0
    CHM.CSi = CSio + (NUM.theta.*dCSidt + (1-NUM.theta).*dCSidto).*NUM.dt;
    CHM.CFe = CFeo + (NUM.theta.*dCFedt + (1-NUM.theta).*dCFedto).*NUM.dt;
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
    CHM.XSi         = MAT.rhot - CHM.XFe;
    CHM.xFe         = CHM.XFe./MAT.rhot;    CHM.xSi     = 1-CHM.xFe;
    % convert chemistry to concentrations (equation 10)
    CHM.cSi         = CHM.CSi./CHM.xSi./MAT.rhot;
    CHM.cFe         = CHM.CFe./CHM.xFe./MAT.rhot;
    % enthalpy/equilibrium method update chemistry
    SOL.T = alpha.*SOL.T + (1-alpha).*SOL.H./(MAT.rhoCpt + MAT.rhot.*CHM.xSi.*CHM.fSil.*CHM.dEntrSi + MAT.rhot.*CHM.xFe.*CHM.fFel.*CHM.dEntrFe);

else
    rhoCpt  = (MAT.rhot.*CHM.xFe.*(CHM.fFes.*PHY.CpFes + CHM.fFel.*PHY.CpFel)...
            +  MAT.rhot.*CHM.xSi.*(CHM.fSis.*PHY.CpSis + CHM.fSil.*PHY.CpSil));
    SOL.H   =  SOL.T  .*(MAT.rhot.*CHM.xFe.*(CHM.fFel.*CHM.dEntrSi) + MAT.rhot.*CHM.xSi.*(CHM.fSil.*CHM.dEntrSi) +rhoCpt);
    CHM.XFe =  MAT.rhot.*CHM.xFe; CHM.XSi = MAT.rhot.*CHM.xSi;
    CHM.CFe =  MAT.rhot.*CHM.cFe.*CHM.xFe;   
    CHM.CSi =  MAT.rhot.*CHM.cSi.*CHM.xSi;  
end

% [CHM.fFes,CHM.csFe,CHM.clFe]     = equilibrium((SOL.T+To)./2,(CHM.cFe+cFeo)./2,(SOL.Pt+Pto)./2,CHM.TFe1,CHM.TFe2,CHM.cphsFe1,CHM.cphsFe2,...
%     CHM.perTFe,CHM.perCsFe,CHM.perClFe,CHM.clap,CHM.PhDgFe,TINY);
%
% [CHM.fSis,CHM.csSi,CHM.clSi]     = equilibrium((SOL.T+To)./2,(CHM.cSi+cSio)./2,(SOL.Pt+Pto)./2,CHM.TSi1,CHM.TSi2,CHM.cphsSi1,CHM.cphsSi2,...
%     CHM.perTSi,CHM.perCsSi,CHM.perClSi,CHM.clap,CHM.PhDgSi,TINY);
% CHM.fFel = 1-CHM.fFes; CHM.fSil = 1-CHM.fSis;

[fFeso,CHM.csFe,CHM.clFe]     = equilibrium(SOL.T,CHM.cFe,SOL.Pt,CHM.TFe1,CHM.TFe2,CHM.cphsFe1,CHM.cphsFe2,...
    CHM.perTFe,CHM.perCsFe,CHM.perClFe,CHM.clap,CHM.PhDgFe,TINY);

[fSiso,CHM.csSi,CHM.clSi]     = equilibrium(SOL.T,CHM.cSi,SOL.Pt,CHM.TSi1,CHM.TSi2,CHM.cphsSi1,CHM.cphsSi2,...
    CHM.perTSi,CHM.perCsSi,CHM.perClSi,CHM.clap,CHM.PhDgSi,TINY);
% lag the phase outputs
CHM.fFes = alpha.*CHM.fFes + (1-alpha).*fFeso; 
% lag the phase outputs
CHM.fSis = alpha.*CHM.fSis + (1-alpha).*fSiso; 

CHM.fFel = 1-CHM.fFes; CHM.fSil = 1-CHM.fSis;

pause(0.1)


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
% [advn]              = advection(SOL.phi,SOL.U,WFe,NUM.h,NUM.h,NUM.AdvnScheme,'Flx');
% if NUM.step>0
% SOL.phi             = phio - advn * NUM.dt;
% end
% SOL.phi(SOL.phi<0)  = 0;
% %% Solve energy equation explicitly
%
% % get advection rate
% advn_T              = advection(SOL.T,SOL.U,SOL.W,NUM.h,NUM.h,NUM.AdvnScheme,'adv');
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