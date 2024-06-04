% planetesimal: solve thermo-chemical equations

if step>0
    step;
end


% test smaller tiny criterion
TINY1 = 1e-6;

% store previous iteration (diagnostic)
Ti    = T;
Si    = S;
xFei  = xFe;
xSii  = xSi;
cFei  = cFe;
cSii  = cSi;
flFei = flFe;
flSii = flSi;

%% update heat content (entropy)
advn_S = - advect(FsSi(inz,inx).*ssSi(inz,inx),UsSi(inz,:),WsSi(:,inx),h,{ADVN,''},[1,2],BCA) ...  % heat advection
         - advect(FsFe(inz,inx).*ssFe(inz,inx),UsFe(inz,:),WsFe(:,inx),h,{ADVN,''},[1,2],BCA) ...
         - advect(FlSi(inz,inx).*slSi(inz,inx),UlSi(inz,:),WlSi(:,inx),h,{ADVN,''},[1,2],BCA) ...
         - advect(FlFe(inz,inx).*slFe(inz,inx),UlFe(inz,:),WlFe(:,inx),h,{ADVN,''},[1,2],BCA);

qSz    = - (ks(1:end-1,:)+ks(2:end,:))./2 .* ddz(T,h);
qSx    = - (ks(:,1:end-1)+ks(:,2:end))./2 .* ddx(T,h);
diff_S = - ddz(qSz(:,2:end-1),h)  ...                                 % heat diffusion
         - ddx(qSx(2:end-1,:),h);

diss_T = EntProd./T(2:end-1,2:end-1);

dSdt   = advn_S + diff_S + diss_T;

% residual of entropy evolution
res_S = (a1*S(inz,inx)-a2*So(inz,inx)-a3*Soo(inz,inx))/dt - (b1*dSdt + b2*dSdto + b3*dSdtoo);

% update solution
S(inz,inx)     = S(inz,inx) - alpha*res_S*dt/a1 + beta*upd_S;
upd_S =   - alpha*res_S*dt/a1 + beta*upd_S;

%% test whether forcing Si to stay constant 

% apply boundaries
Ds = xFe.*fsFe.*dEntrFe + xSi.*fsSi.*dEntrSi;

switch BCTTop
    case 'isothermal'
        S(1,:)  = RHO(1,:).*(Cp.*log(Ttop0./T0)+aT./rhoRef.*(Pt(1,:)-P0) + Ds(1,:));
    case 'insulating'
        S(1,:)  = S(2,:);
    case 'flux'
        T2      = T0.*exp((S(2,:) - FsFe(2,:).*dEntrFe - FsSi(2,:).*dEntrSi)./RHO(2,:)./Cp ...
                + aT.*(Pt(2,:) - P0)./rhoRef./Cp);
        dTb     = qT0*h./((kT(1,:)+kT(2,:))./2);
        Tb      = max(Ttop0,T2+dTb);
        S(1,:)  = RHO(1,:).*(Cp.*log(Tb./T0)+aT./rhoRef.*(Pt(1,:)-P0) + Ds(1,:));
end
switch BCTBot
    case 'isothermal'
S(end,:)  = RHO(end,:).*(Cp.*log(Tbot0./T0)+aT./rhoRef.*(Pt(end,:)-P0) + Ds(end,:));    
    case 'insulating'
        S(end,:) = S(end-1,:);
end
switch BCTSides
    case 'isothermal'
        S(:,[1 end]) = RHO(:,[1 end]).*(Cp.*log(T0./T0)+aT./rhoRef.*(Pt(:,[1 end])-P0) + Ds(:,[1 end]));
    case 'insulating'
        S(:,[1 end]) = S(:,[2 end-1]);
end

%% update system fractions, only one system needs to be solved as SUM_i(X_i) = 1
if any(xFe(:)>0 & xFe(:)<1)
    dXFedt          = - advect(FsFe(inz,inx),UsFe(inz,:),WsFe(:,inx),h,{ADVN,''},[1,2],BCA) ...
                      - advect(FlFe(inz,inx),UlFe(inz,:),WlFe(:,inx),h,{ADVN,''},[1,2],BCA);

    % update solution
    XFe(inz,inx)    = (a2*XFeo(inz,inx) + a3*XFeoo(inz,inx) + (b1*dXFedt + b2*dXFedto + b3*dXFedtoo)*dt)/a1;


    dXSidt          = - advect(FsSi(inz,inx),UsSi(inz,:),WsSi(:,inx),h,{ADVN,''},[1,2],BCA) ...
                      - advect(FlSi(inz,inx),UlSi(inz,:),WlSi(:,inx),h,{ADVN,''},[1,2],BCA);
    

    % update solution
    % residual of entropy evolution
    res_XFe = (a1*XFe(inz,inx)-a2*XFeo(inz,inx)-a3*XFeoo(inz,inx))/dt - (b1*dXFedt + b2*dXFedto + b3*dXFedtoo);
    res_XSi = (a1*XSi(inz,inx)-a2*XSio(inz,inx)-a3*XSioo(inz,inx))/dt - (b1*dXSidt + b2*dXSidto + b3*dXSidtoo);

    % update solution
    % semi-implicit update of bulk entropy density
    XFe(inz,inx)     = XFe(inz,inx) - alpha*res_XFe*dt/a1 + beta*upd_XFe;
    XSi(inz,inx)     = XSi(inz,inx) - alpha*res_XSi*dt/a1 + beta*upd_XSi;
    upd_XFe =   - alpha*res_XFe*dt/a1 + beta*upd_XFe;
    upd_XSi =   - alpha*res_XSi*dt/a1 + beta*upd_XSi;

    % apply boundaries
    XFe([1 end],:)  = XFe([2 end-1],:);  XFe(:,[1 end]) = XFe(:,[2 end-1]);
    XSi([1 end],:)  = XSi([2 end-1],:);  XSi(:,[1 end]) = XSi(:,[2 end-1]);

    % enforce 0,rho limits
    XFe = max(0, XFe ); 
    XSi = max(0, XSi ); 

else
    XFe = xFe.*rho;
    XSi = rho - XFe;
end
    RHO = XFe+XSi;

%% update composition

% update fertile chemical components
% equations 8a, 8b

advn_CSi = - advect(FsSi(inz,inx).*csSi(inz,inx),UsSi(inz,:),WsSi(:,inx),h,{ADVN,''},[1,2],BCA) ...
           - advect(FlSi(inz,inx).*clSi(inz,inx),UlSi(inz,:),WlSi(:,inx),h,{ADVN,''},[1,2],BCA);

dCSidt   = advn_CSi;

advn_CFe = - advect(FsFe(inz,inx).*csFe(inz,inx),UsFe(inz,:),WsFe(:,inx),h,{ADVN,''},[1,2],BCA) ...
           - advect(FlFe(inz,inx).*clFe(inz,inx),UlFe(inz,:),WlFe(:,inx),h,{ADVN,''},[1,2],BCA);

dCFedt   = advn_CFe;

% update solution
% residual of entropy evolution
res_CFe = (a1*CFe(inz,inx)-a2*CFeo(inz,inx)-a3*CFeoo(inz,inx))/dt - (b1*dCFedt + b2*dCFedto + b3*dCFedtoo);
res_CSi = (a1*CSi(inz,inx)-a2*CSio(inz,inx)-a3*CSioo(inz,inx))/dt - (b1*dCSidt + b2*dCSidto + b3*dCSidtoo);

% update solution
% semi-implicit update of bulk chemical composition density
CFe(inz,inx)     = CFe(inz,inx) - alpha*res_CFe*dt/a1 + beta*upd_CFe;
CSi(inz,inx)     = CSi(inz,inx) - alpha*res_CSi*dt/a1 + beta*upd_CSi;
upd_CFe =   - alpha*res_CFe*dt/a1 + beta*upd_CFe;
upd_CSi =   - alpha*res_CSi*dt/a1 + beta*upd_CSi;

% cheat a little bit again and force C_i = X_i*c_i(t-1) when below the
% solidus or above the liquidus
CSi(~hassolSi|~hasliqSi) = XSi(~hassolSi|~hasliqSi) .*cSi(~hassolSi|~hasliqSi);
CFe(~hassolFe|~hasliqFe) = XFe(~hassolFe|~hasliqFe) .*cFe(~hassolFe|~hasliqFe);

% apply boundaries
CSi([1 end],:) = CSi([2 end-1],:);  CSi(:,[1 end]) = CSi(:,[2 end-1]);
CFe([1 end],:) = CFe([2 end-1],:);  CFe(:,[1 end]) = CFe(:,[2 end-1]);

% % enforce 0,rho limits
CFe = max(0,CFe);
CSi = max(0,CSi);

%% update local phase equilibrium
[fsFeq,csFeq,clFeq] = equilibrium(T,cFe,Pt,TFe1,TFe2,cphsFe1,cphsFe2,...
                                  perTFe,percsFe,perclFe,clap,PhDgFe);
% apply boundaries
fsFeq([1 end],:) = fsFeq([2 end-1],:);
fsFeq(:,[1 end]) = fsFeq(:,[2 end-1]);
csFeq([1 end],:) = csFeq([2 end-1],:);
csFeq(:,[1 end]) = csFeq(:,[2 end-1]);
clFeq([1 end],:) = clFeq([2 end-1],:);
clFeq(:,[1 end]) = clFeq(:,[2 end-1]);

[fsSiq,csSiq,clSiq] = equilibrium(T,cSi,Pt,TSi1,TSi2,cphsSi1,cphsSi2,...
                                  perTSi,percsSi,perclSi,clap,PhDgSi);
% apply boundaries
fsSiq([1 end],:) = fsSiq([2 end-1],:);
fsSiq(:,[1 end]) = fsSiq(:,[2 end-1]);
csSiq([1 end],:) = csSiq([2 end-1],:);
csSiq(:,[1 end]) = csSiq(:,[2 end-1]);
clSiq([1 end],:) = clSiq([2 end-1],:);
clSiq(:,[1 end]) = clSiq(:,[2 end-1]);

flFeq = 1-fsFeq;
flSiq = 1-fsSiq;

%% update phase fractions
% solid
GFes        = ((XFe.*fsFeq-FsFe)./(4.*dt));
advn_FFes   = - advect(FsFe(inz,inx),UsFe(inz,:),WsFe(:,inx),h,{ADVN,''},[1,2],BCA);
dFsFedt     = advn_FFes + GFes(inz,inx);
res_FsFe = (a1*FsFe(inz,inx)-a2*FsFeo(inz,inx)-a3*FsFeoo(inz,inx))/dt - (b1*dFsFedt + b2*dFsFedto + b3*dFsFedtoo);

GSis        = ((XSi.*fsSiq-FsSi)./(4.*dt));
advn_FSis   = - advect(FsSi(inz,inx),UsSi(inz,:),WsSi(:,inx),h,{ADVN,''},[1,2],BCA);
dFsSidt     = advn_FSis + GSis(inz,inx);                                       % total rate of change
res_FsSi = (a1*FsSi(inz,inx)-a2*FsSio(inz,inx)-a3*FsSioo(inz,inx))/dt - (b1*dFsSidt + b2*dFsSidto + b3*dFsSidtoo);

FsFe(inz,inx)     = FsFe(inz,inx) - alpha*res_FsFe*dt/a1 + beta*upd_FsFe;
FsSi(inz,inx)     = FsSi(inz,inx) - alpha*res_FsSi*dt/a1 + beta*upd_FsSi;
upd_FsFe =   - alpha*res_FsFe*dt/a1 + beta*upd_FsFe;
upd_FsSi =   - alpha*res_FsSi*dt/a1 + beta*upd_FsSi;

FsFe([1 end],:) = FsFe([2 end-1],:);  FsFe(:,[1 end]) = FsFe(:,[2 end-1]);  % apply boundary conditions
FsSi([1 end],:) = FsSi([2 end-1],:);  FsSi(:,[1 end]) = FsSi(:,[2 end-1]);  % apply boundary conditions

FsFe = max(0, FsFe );
FsSi = max(0, FsSi );

% liquid
GFel        = ((XFe.*flFeq-FlFe)./(4.*dt));
% GFel = -GFes;
advn_FFel   = - advect(FlFe(inz,inx),UlFe(inz,:),WlFe(:,inx),h,{ADVN,''},[1,2],BCA);
dFlFedt     = advn_FFel + GFel(inz,inx);
res_FlFe = (a1*FlFe(inz,inx)-a2*FlFeo(inz,inx)-a3*FlFeoo(inz,inx))/dt - (b1*dFlFedt + b2*dFlFedto + b3*dFlFedtoo);

GSil        = ((XSi.*flSiq-FlSi)./(4.*dt));
% GSil = -GSis;
advn_FSil   = - advect(FlSi(inz,inx),UlSi(inz,:),WlSi(:,inx),h,{ADVN,''},[1,2],BCA);
dFlSidt     = advn_FSil + GSil(inz,inx);                                       % total rate of change
res_FlSi = (a1*FlSi(inz,inx)-a2*FlSio(inz,inx)-a3*FlSioo(inz,inx))/dt - (b1*dFlSidt + b2*dFlSidto + b3*dFlSidtoo);

FlFe(inz,inx)     = FlFe(inz,inx) - alpha*res_FlFe*dt/a1 + beta*upd_FlFe;
FlSi(inz,inx)     = FlSi(inz,inx) - alpha*res_FlSi*dt/a1 + beta*upd_FlSi;

upd_FlFe =   - alpha*res_FlFe*dt/a1 + beta*upd_FlFe;
upd_FlSi =   - alpha*res_FlSi*dt/a1 + beta*upd_FlSi;
FlFe([1 end],:) = FlFe([2 end-1],:);  FlFe(:,[1 end]) = FlFe(:,[2 end-1]);  % apply boundary conditions
FlSi([1 end],:) = FlSi([2 end-1],:);  FlSi(:,[1 end]) = FlSi(:,[2 end-1]);  % apply boundary conditions

FlFe = max(0, FlFe );
FlSi = max(0, FlSi );

RHO = FsSi+FsFe+FlSi+FlFe;

if radheat 
% update radioactive isotope decay and heating rate
[dndt,H]    = rad_decay(n26Al,tauAl,EAl);
Alfrac      = (cSi-cphsSi1)/(cphsSi2-cphsSi1);
Hr          = H.*XSi.*Alfrac;
n26Al = (a2*n26Alo + a3*n26Aloo + (b1*dndt + b2*dndto + b3*dndtoo)*dt)/a1;
end

% % update temperature
T   = T0.*exp((S - FsFe.*dEntrFe - FsSi.*dEntrSi)./(FlFe+FsFe+FlSi+FsSi)./Cp ...
    + aT.*(Pt - P0)./rhoRef./Cp);

% update system fractions
xFe = max(0,min(1, XFe./RHO));
xSi = max(0,min(1, XSi./RHO ));

hasFe   = xFe>TINY1 & xSi<1-TINY1;
hasSi   = xSi>TINY1 & xFe<1-TINY1;

% update chemical composition
cSi(hasSi) = CSi(hasSi)./max(XSi(hasSi),(FlSi(hasSi)+FsSi(hasSi)));
cFe(hasFe) = CFe(hasFe)./max(XFe(hasFe),(FlFe(hasFe)+FsFe(hasFe)));

% update phase fractions [wt]
flFe(hasFe) = max(0,min(1,FlFe(hasFe)./XFe(hasFe)));
flSi(hasSi) = max(0,min(1,FlSi(hasSi)./XSi(hasSi)));
fsFe(hasFe) = max(0,min(1,FsFe(hasFe)./XFe(hasFe)));
fsSi(hasSi) = max(0,min(1,FsSi(hasSi)./XSi(hasSi)));


% flFe(hasFe) = 1-fsFe(hasFe); 
% flSi(hasSi) = 1-fsSi(hasSi);

%detect where is fully molten
hassolSi = flSi<1-TINY1 & fsSi>TINY1;
hassolFe = flFe<1-TINY1 & fsFe>TINY1;

hasliqSi = fsSi<1-TINY1 & flSi>TINY1;
hasliqFe = fsFe<1-TINY1 & flFe>TINY1;

%% update phase compositions
KcFe = csFeq./clFeq;
clFe(hasFe) = cFe(hasFe)./max(TINY,(flFe(hasFe) + fsFe(hasFe).*KcFe(hasFe)));
csFe(hasFe) = cFe(hasFe)./max(TINY,(flFe(hasFe)./KcFe(hasFe) + fsFe(hasFe)));

KcSi = csSiq./clSiq;
clSi(hasSi) = cSi(hasSi)./max(TINY,(flSi(hasSi) + fsSi(hasSi).*KcSi(hasSi)));
csSi(hasSi) = cSi(hasSi)./max(TINY,(flSi(hasSi)./KcSi(hasSi) + fsSi(hasSi)));

%% update phase entropies
slFe  = (S - FsFe.*dEntrFe - FsSi.*dEntrSi)./(FlFe+FsFe+FlSi+FsSi);
slSi  = slFe;
ssFe  = slFe + dEntrFe;
ssSi  = slSi + dEntrSi;

%% get residual of thermochemical equations from iterative update
normS   = norm(S    - Si   ,2)./(norm(S   ,2)+TINY);
normxFe = norm(xFe  - xFei ,2)./(norm(xFe ,2)+TINY);
normcFe = norm(cFe  - cFei ,2)./(norm(cFe ,2)+TINY);
normcSi = norm(cSi  - cSii ,2)./(norm(cSi ,2)+TINY);
normfFe = norm(flFe - flFei,2)./(norm(flFe,2)+TINY);
normfSi = norm(flSi - flSii,2)./(norm(flSi,2)+TINY);
resnorm_TC = normS + normxFe + normcSi + normcFe + normfFe + normfSi;