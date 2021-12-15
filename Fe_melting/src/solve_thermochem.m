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
advn_H      = advection(XSi.*fSis.*PHY.CpSis.*SOL.T,...
                        UsSi,WsSi,NUM.h,NUM.h,NUM.AdvnScheme,'flx')...
            + advection(XFe.*fFes.*PHY.CpFes.*SOL.T,...
                        UsFe,WsFe,NUM.h,NUM.h,NUM.AdvnScheme,'flx')...
            + advection(XSi.*fSil.*SOL.T.*(PHY.CpSil + dEntrSi),...
                        UlSi,WlSi,NUM.h,NUM.h,NUM.AdvnScheme,'flx')...
            + advection(XFe.*fFel.*SOL.T.*(PHY.CpFel + dEntrFe),...
                        UlFe,WlFe,NUM.h,NUM.h,NUM.AdvnScheme,'flx');

qTz    = - (MAT.kT(1:end-1,:)+MAT.kT(2:end,:))./2 .* ddz(SOL.T,NUM.h);                 % heat diffusion z-flux
qTx    = - (MAT.kT(:,1:end-1)+MAT.kT(:,2:end))./2 .* ddx(SOL.T,NUM.h);                 % heat diffusion x-flux
diff_T(2:end-1,2:end-1) = (- ddz(qTz(:,2:end-1),NUM.h) ...                 % heat diffusion
            - ddx(qTx(2:end-1,:),NUM.h));
dHdt        = diff_T - advn_H +MAT.Hr +zeros(size(SOL.T));

%% update composition

% update system, only one system needs to be solved as SUM_i(X_i) = 1
% equation 7
adv_XFe     = advection(XFe.*fFes,UsFe,WsFe,NUM.h,NUM.h,NUM.AdvnScheme,'flx')...
            + advection(XFe.*fFel,UlFe,WlFe,NUM.h,NUM.h,NUM.AdvnScheme,'flx');
dXdt        = -adv_XFe;

% update fertile chemical components
% equations 8a, 8b
advn_CSi    = advection(XSi.*fSis.*csSi,...
                        UsSi,WsSi,NUM.h,NUM.h,NUM.AdvnScheme,'flx')...
            + advection(XSi.*fSil.*clSi,...
                        UlSi,WlSi,NUM.h,NUM.h,NUM.AdvnScheme,'flx');

qCSiz    = - (MAT.kC(1:end-1,:).*fSil(1:end-1,:).*MAT.rhoSil(1:end-1,:) ...
             +MAT.kC(2:end,:)  .*fSil(2:end,:)  .*MAT.rhoSil(2:end,:))./2 ...
           .* ddz(SOL.T,NUM.h);                 % heat diffusion z-flux
qCSix    = - (MAT.kC(:,1:end-1).*fSil(:,1:end-1).*MAT.rhoSil(:,1:end-1)...
            + MAT.kC(:,2:end)  .*fSil(:,2:end)  .*MAT.rhoSil(:,2:end))./2 ...
           .* ddx(SOL.T,NUM.h);                 % heat diffusion x-flux
diff_CSi(2:end-1,2:end-1) = (- ddz(qCSiz(:,2:end-1),NUM.h) ...                 % heat diffusion
    - ddx(qCSix(2:end-1,:),NUM.h));

dCSidt      = -advn_CSi + diff_CSi +zeros(size(CSi));

advn_CFe    = advection(XFe.*fFes.*csFe,...
                        UsFe,WsFe,NUM.h,NUM.h,NUM.AdvnScheme,'flx')...
            + advection(XFe.*fFel.*clFe,...
                        UlFe,WlFe,NUM.h,NUM.h,NUM.AdvnScheme,'flx');

qCFez    = - (MAT.kC(1:end-1,:).*fFel(1:end-1,:).*MAT.rhoFel(1:end-1,:) ...
             +MAT.kC(2:end,:)  .*fFel(2:end,:)  .*MAT.rhoFel(2:end,:))./2 ...
           .* ddz(SOL.T,NUM.h);                 % heat diffusion z-flux
qCFex    = - (MAT.kC(:,1:end-1).*fFel(:,1:end-1).*MAT.rhoFel(:,1:end-1)...
            + MAT.kC(:,2:end)  .*fFel(:,2:end)  .*MAT.rhoFel(:,2:end))./2 ...
           .* ddx(SOL.T,NUM.h);                 % heat diffusion x-flux
diff_CFe(2:end-1,2:end-1) = (- ddz(qCFez(:,2:end-1),NUM.h) ...                 % heat diffusion
    - ddx(qCFex(2:end-1,:),NUM.h));

dCFedt      = -advn_CFe + diff_CFe + zeros(size(CFe));

if NUM.step>0   
    SOL.H       = Ho + (NUM.theta.*dHdt + (1-NUM.theta).*dHdto).*NUM.dt;
    CSi = CSio + (NUM.theta.*dCSidt + (1-NUM.theta).*dCSidto).*NUM.dt;
    CFe = CFeo + (NUM.theta.*dCFedt + (1-NUM.theta).*dCFedto).*NUM.dt;
    XFe = XFeo + (NUM.theta.*dXdt   + (1-NUM.theta).*dXdto)  .*NUM.dt;

XSi         = MAT.rhot - XFe;
xFe         = XFe./MAT.rhot;    xSi     = 1-xFe;
% convert chemistry to concentrations (equation 9)
cSi         = CSi./XSi./MAT.rhot;
cFe         = CFe./XFe./MAT.rhot;
end
% enthalpy/equilibrium method update chemistry
SOL.T = SOL.H./(MAT.rhoCpt + XSi.*fSil.*dEntrSi + XFe.*fFel.*dEntrFe);


[fFes,csFe,clFe]     = equilibrium((SOL.T+To)./2,(cFe+cFeo)./2,(SOL.Pt+Pto)./2,TFe1,TFe2,cphsFe1,cphsFe2,...
    perTSi,perCsFe,perClFe,clap,PhDgFe,TINY);
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