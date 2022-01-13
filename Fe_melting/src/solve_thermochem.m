% planetesimal: solve thermo-chemical equations

% print solver header
fprintf(1,'  ---  solve thermo-chemical equations \n');
tic;

bndmode = 0;

switch bndmode
    case 0
        bndshape = zeros(size(SOL.T));
    case 1
        bndshape = exp( ( -ZZ)/dw);
    case 2
        bndshape = exp( ( -ZZ)/dw) ...
                 + exp(-(D-ZZ)/dw);
    case 3
        bndshape = exp( ( -ZZ)/dw) ...
                 + exp(-(D-ZZ)/dw) ...
                 + exp( ( -XX)/dw) ...
                 + exp(-(L-XX)/dw);
end


%% update enthalpy and Temperature
advn_H      = advection(MAT.rhot.*xSi.*fSis.*PHY.CpSis.*SOL.T,...
                        UsSi,WsSi,NUM.h,NUM.h,NUM.AdvnScheme,'flx')...
            + advection(MAT.rhot.*xFe.*fFes.*PHY.CpFes.*SOL.T,...
                        UsFe,WsFe,NUM.h,NUM.h,NUM.AdvnScheme,'flx')...
            + advection(MAT.rhot.*xSi.*fSil.*SOL.T.*(PHY.CpSil + dEntrSi),...
                        UlSi,WlSi,NUM.h,NUM.h,NUM.AdvnScheme,'flx')...
            + advection(MAT.rhot.*xFe.*fFel.*SOL.T.*(PHY.CpFel + dEntrFe),...
                        UlFe,WlFe,NUM.h,NUM.h,NUM.AdvnScheme,'flx');

qTz    = - (MAT.kT(1:end-1,:)+MAT.kT(2:end,:))./2 .* ddz(SOL.T,NUM.h);                 % heat diffusion z-flux
qTx    = - (MAT.kT(:,1:end-1)+MAT.kT(:,2:end))./2 .* ddx(SOL.T,NUM.h);                 % heat diffusion x-flux
diff_T(2:end-1,2:end-1) = (- ddz(qTz(:,2:end-1),NUM.h) ...                 % heat diffusion
            - ddx(qTx(2:end-1,:),NUM.h));
dHdt        = diff_T - advn_H +MAT.Hr +zeros(size(SOL.T));

%% update composition

% update system, only one system needs to be solved as SUM_i(X_i) = 1
% equation 7
adv_XFe     = advection(MAT.rhot.*xFe.*fFes,UsFe,WsFe,NUM.h,NUM.h,NUM.AdvnScheme,'flx')...
            + advection(MAT.rhot.*xFe.*fFel,UlFe,WlFe,NUM.h,NUM.h,NUM.AdvnScheme,'flx');
dXdt        = -adv_XFe;

% update fertile chemical components
% equations 8a, 8b
advn_CSi    = advection(MAT.rhot.*xSi.*fSis.*csSi,...
                        UsSi,WsSi,NUM.h,NUM.h,NUM.AdvnScheme,'flx')...
            + advection(MAT.rhot.*xSi.*fSil.*clSi,...
                        UlSi,WlSi,NUM.h,NUM.h,NUM.AdvnScheme,'flx');

qCSiz    = - (MAT.kC(1:end-1,:).*fSil(1:end-1,:).*MAT.rhot(1:end-1,:) ...
             +MAT.kC(2:end,:)  .*fSil(2:end,:)  .*MAT.rhot(2:end,:))./2 ...
           .* ddz(clSi,NUM.h);                 % heat diffusion z-flux
qCSix    = - (MAT.kC(:,1:end-1).*fSil(:,1:end-1).*MAT.rhot(:,1:end-1)...
            + MAT.kC(:,2:end)  .*fSil(:,2:end)  .*MAT.rhot(:,2:end))./2 ...
           .* ddx(clSi,NUM.h);                 % heat diffusion x-flux
diff_CSi(2:end-1,2:end-1) = (- ddz(qCSiz(:,2:end-1),NUM.h) ...                 % heat diffusion
    - ddx(qCSix(2:end-1,:),NUM.h));

dCSidt      = -advn_CSi + diff_CSi +zeros(size(CSi));

advn_CFe    = advection(MAT.rhot.*xFe.*fFes.*csFe,...
                        UsFe,WsFe,NUM.h,NUM.h,NUM.AdvnScheme,'flx')...
            + advection(MAT.rhot.*xFe.*fFel.*clFe,...
                        UlFe,WlFe,NUM.h,NUM.h,NUM.AdvnScheme,'flx');

qCFez    = - (MAT.kC(1:end-1,:).*fFel(1:end-1,:).*MAT.rhot(1:end-1,:) ...
             +MAT.kC(2:end,:)  .*fFel(2:end,:)  .*MAT.rhot(2:end,:))./2 ...
           .* ddz(clFe,NUM.h);                 % heat diffusion z-flux
qCFex    = - (MAT.kC(:,1:end-1).*fFel(:,1:end-1).*MAT.rhot(:,1:end-1)...
            + MAT.kC(:,2:end)  .*fFel(:,2:end)  .*MAT.rhot(:,2:end))./2 ...
           .* ddx(clFe,NUM.h);                 % heat diffusion x-flux
diff_CFe(2:end-1,2:end-1) = (- ddz(qCFez(:,2:end-1),NUM.h) ...                 % heat diffusion
    - ddx(qCFex(2:end-1,:),NUM.h));

dCFedt      = -advn_CFe + diff_CFe + zeros(size(CFe));

if NUM.step>0   
    SOL.H = Ho + (NUM.theta.*dHdt   + (1-NUM.theta).*dHdto)  .*NUM.dt;    SOL.H([1 end],:) = SOL.H([2 end-1],:);  SOL.H(:,[1 end]) = SOL.H(:,[2 end-1]);
    CSi = CSio + (NUM.theta.*dCSidt + (1-NUM.theta).*dCSidto).*NUM.dt;    CSi([1 end],:) = CSi([2 end-1],:);  CSi(:,[1 end]) = CSi(:,[2 end-1]);
    CFe = CFeo + (NUM.theta.*dCFedt + (1-NUM.theta).*dCFedto).*NUM.dt;    CFe([1 end],:) = CFe([2 end-1],:);  CFe(:,[1 end]) = CFe(:,[2 end-1]);
    XFe = XFeo + (NUM.theta.*dXdt   + (1-NUM.theta).*dXdto)  .*NUM.dt;    XFe([1 end],:) = XFe([2 end-1],:);  XFe(:,[1 end]) = XFe(:,[2 end-1]);

XSi         = MAT.rhot - XFe;
xFe         = XFe./MAT.rhot;    xSi     = 1-xFe;
% convert chemistry to concentrations (equation 9)
cSi         = alpha.*cSi + (1-alpha).*CSi./xSi./MAT.rhot;
cFe         = alpha.*cFe + (1-alpha).*CFe./xFe./MAT.rhot;
end
% enthalpy/equilibrium method update chemistry
SOL.T = alpha.*SOL.T + (1-alpha).*SOL.H./(MAT.rhoCpt + MAT.rhot.*xSi.*fSil.*dEntrSi + MAT.rhot.*xFe.*fFel.*dEntrFe);


[fFes,csFe,clFe]     = equilibrium((SOL.T+To)./2,(cFe+cFeo)./2,(SOL.Pt+Pto)./2,TFe1,TFe2,cphsFe1,cphsFe2,...
    perTFe,perCsFe,perClFe,clap,PhDgFe,TINY);

[fSis,csSi,clSi]     = equilibrium((SOL.T+To)./2,(cSi+cSio)./2,(SOL.Pt+Pto)./2,TSi1,TSi2,cphsSi1,cphsSi2,...
    perTSi,perCsSi,perClSi,clap,PhDgSi,TINY);
fFel = 1-fFes; fSil = 1-fSis;


% %check outputs
% NUM.dt
% 
% figure(100)
% subplot(2,2,1)
% imagesc(SOL.H-Ho); colorbar
% title('dH')
% 
% subplot(2,2,2)
% imagesc(XFe - XFeo); colorbar
% title('dX')
% 
% subplot(2,2,3)
% imagesc(CSi - CSio); colorbar
% title('dCSi')
% 
% subplot(2,2,4)
% imagesc(CFe - CFeo); colorbar
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