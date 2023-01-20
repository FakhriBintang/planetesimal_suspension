% planetesimal: solve thermo-chemical equations

if NUM.step>0
    NUM.step;
end

% store previous iteration
Ti    = SOL.T;
xFei  = CHM.xFe;
xSii  = CHM.xSi;
cFei  = CHM.cFe;
cSii  = CHM.cSi;
flFei = CHM.flFe;
flSii = CHM.flSi;


%% update phase entropies
SOL.slFe  = (SOL.S - CHM.FsFe.*CHM.dEntrFe - CHM.FsSi.*CHM.dEntrSi)./MAT.rho;
SOL.slSi  = SOL.slFe;
SOL.ssFe  = SOL.slFe + CHM.dEntrFe;
SOL.ssSi  = SOL.slSi + CHM.dEntrSi;


%% update heat content (entropy)
advn_S = - advect(CHM.FsSi(inz,inx).*SOL.ssSi(inz,inx),UsSi(inz,:),WsSi(:,inx),NUM.h,{ADVN,''},[1,2],BCA) ...  % heat advection
         - advect(CHM.FsFe(inz,inx).*SOL.ssFe(inz,inx),UsFe(inz,:),WsFe(:,inx),NUM.h,{ADVN,''},[1,2],BCA) ...
         - advect(CHM.FlSi(inz,inx).*SOL.slSi(inz,inx),UlSi(inz,:),WlSi(:,inx),NUM.h,{ADVN,''},[1,2],BCA) ...
         - advect(CHM.FlFe(inz,inx).*SOL.slFe(inz,inx),UlFe(inz,:),WlFe(:,inx),NUM.h,{ADVN,''},[1,2],BCA);

qSz    = - (MAT.ks(1:end-1,:)+MAT.ks(2:end,:))./2 .* ddz(SOL.T,NUM.h);
qSx    = - (MAT.ks(:,1:end-1)+MAT.ks(:,2:end))./2 .* ddx(SOL.T,NUM.h);
diff_S = - ddz(qSz(:,2:end-1),NUM.h)  ...                                 % heat diffusion
         - ddx(qSx(2:end-1,:),NUM.h);

diss_T = EntProd./SOL.T(2:end-1,2:end-1);

dSdt   = advn_S + diff_S + MAT.Hr + diss_T;

% update solution
SOL.S(inz,inx) = So(inz,inx) + (NUM.theta.*dSdt + (1-NUM.theta).*dSdto) .* NUM.dt;

% apply boundaries
MAT.Ds = CHM.xFe.*CHM.fsFe.*CHM.dEntrFe + CHM.xSi.*CHM.fsSi.*CHM.dEntrSi;
switch SOL.BCTTop
    case 'isothermal'
        SOL.S(1,:) = MAT.rho(1,:).*(PHY.Cp.*log(SOL.T0./SOL.T0)+PHY.aT./rhoRef.*(SOL.Pt(1,:)-P0) + MAT.Ds(1,:));
    case 'insulating'
        SOL.S(1,:) = SOL.S(2,:);
    case 'flux'
end
switch SOL.BCTBot
    case 'isothermal'
        SOL.S(end,:) = MAT.rho(end,:).*(PHY.Cp.*log(SOL.T1./SOL.T0)+PHY.aT./rhoRef.*(SOL.Pt(end,:)-P0) + MAT.Ds(end,:));
    case 'insulating'
        SOL.S(end,:) = SOL.S(end-1,:);
end
switch SOL.BCTSides
    case 'isothermal'
        SOL.S(:,[1 end]) = MAT.rho(:,[1 end]).*(PHY.Cp.*log(SOL.T0./SOL.T0)+PHY.aT./rhoRef.*(SOL.Pt(:,[1 end])-P0) + MAT.Ds(:,[1 end]));
    case 'insulating'
        SOL.S(:,[1 end]) = SOL.S(:,[2 end-1]);
end

% update temperature
SOL.T   = SOL.T0.*exp((SOL.S - CHM.FsFe.*CHM.dEntrFe - CHM.FsSi.*CHM.dEntrSi)./MAT.rho./PHY.Cp ...
        + PHY.aT.*(SOL.Pt - P0)./rhoRef./PHY.Cp);


%% update system fractions, only one system needs to be solved as SUM_i(X_i) = 1
% equation 7

% update system indices
hasFe   = CHM.xFe         >0;
hasSi   = CHM.xSi         >0;

% if any(CHM.xFe(:)>0 & CHM.xFe(:)<1)
    dXFedt = - advect(CHM.FsFe(inz,inx),UsFe(inz,:),WsFe(:,inx),NUM.h,{ADVN,''},[1,2],BCA) ...
             - advect(CHM.FlFe(inz,inx),UlFe(inz,:),WlFe(:,inx),NUM.h,{ADVN,''},[1,2],BCA);

    % update solution
    CHM.XFe(inz,inx) = XFeo(inz,inx) + (NUM.theta.*dXFedt + (1-NUM.theta).*dXFedto) .* NUM.dt;

    dXSidt = - advect(CHM.FsSi(inz,inx),UsSi(inz,:),WsSi(:,inx),NUM.h,{ADVN,''},[1,2],BCA) ...
             - advect(CHM.FlSi(inz,inx),UlSi(inz,:),WlSi(:,inx),NUM.h,{ADVN,''},[1,2],BCA);
    
    % update solution
    CHM.XSi(inz,inx) = XSio(inz,inx) + (NUM.theta.*dXSidt + (1-NUM.theta).*dXSidto) .* NUM.dt;

    CHM.XFe(~hasFe) = MAT.rho(~hasFe) - CHM.XSi(~hasFe);
    CHM.XSi( hasFe) = MAT.rho( hasFe) - CHM.XFe( hasFe);

    % apply boundaries
    CHM.XFe([1 end],:) = CHM.XFe([2 end-1],:);  CHM.XFe(:,[1 end]) = CHM.XFe(:,[2 end-1]);
    CHM.XSi([1 end],:) = CHM.XSi([2 end-1],:);  CHM.XSi(:,[1 end]) = CHM.XSi(:,[2 end-1]);

    % enforce 0,rho limits
    CHM.XFe = max(0, CHM.XFe ); 
    CHM.XSi = max(0, CHM.XSi ); 

% else
%     CHM.XFe = CHM.xFe.*MAT.rho;
%     CHM.XSi = MAT.rho - CHM.XFe;
% end

% update system fractions
CHM.xFe = max(0,min(1, CHM.XFe./MAT.rho ));
CHM.xSi = max(0,min(1, CHM.XSi./MAT.rho ));

hasFe   = CHM.xFe         >0;
hasSi   = CHM.xSi         >0;
hasFein = CHM.xFe(inz,inx)>0;
hasSiin = CHM.xSi(inz,inx)>0;


%% update composition

% update fertile chemical components
% equations 8a, 8b

advn_CSi = - advect(CHM.FsSi(inz,inx).*CHM.csSi(inz,inx),UsSi(inz,:),WsSi(:,inx),NUM.h,{ADVN,''},[1,2],BCA) ...
           - advect(CHM.FlSi(inz,inx).*CHM.clSi(inz,inx),UlSi(inz,:),WlSi(:,inx),NUM.h,{ADVN,''},[1,2],BCA);

dCSidt   = advn_CSi;

advn_CFe = - advect(CHM.FsFe(inz,inx).*CHM.csFe(inz,inx),UsFe(inz,:),WsFe(:,inx),NUM.h,{ADVN,''},[1,2],BCA) ...
           - advect(CHM.FlFe(inz,inx).*CHM.clFe(inz,inx),UlFe(inz,:),WlFe(:,inx),NUM.h,{ADVN,''},[1,2],BCA);

dCFedt   = advn_CFe;

% update solution
CHM.CSi(inz,inx) = CSio(inz,inx) + (NUM.theta.*dCSidt + (1-NUM.theta).*dCSidto) .* NUM.dt;
CHM.CFe(inz,inx) = CFeo(inz,inx) + (NUM.theta.*dCFedt + (1-NUM.theta).*dCFedto) .* NUM.dt;

% apply boundaries
CHM.CSi([1 end],:) = CHM.CSi([2 end-1],:);  CHM.CSi(:,[1 end]) = CHM.CSi(:,[2 end-1]);
CHM.CFe([1 end],:) = CHM.CFe([2 end-1],:);  CHM.CFe(:,[1 end]) = CHM.CFe(:,[2 end-1]);

% enforce 0,rho limits
CHM.CFe = max(0,CHM.CFe);
CHM.CSi = max(0,CHM.CSi);

% update chemical composition
CHM.cSi(hasSi) = CHM.CSi(hasSi)./CHM.XSi(hasSi);
CHM.cFe(hasFe) = CHM.CFe(hasFe)./CHM.XFe(hasFe);


%% update local phase equilibrium
[fsFeq,csFeq,clFeq] = equilibrium(SOL.T,CHM.cFe,SOL.Pt,CHM.TFe1,CHM.TFe2,CHM.cphsFe1,CHM.cphsFe2,...
                                  CHM.perTFe,CHM.perCsFe,CHM.perClFe,CHM.clap,CHM.PhDgFe);
% apply boundaries
fsFeq([1 end],:) = fsFeq([2 end-1],:);
fsFeq(:,[1 end]) = fsFeq(:,[2 end-1]);
csFeq([1 end],:) = csFeq([2 end-1],:);
csFeq(:,[1 end]) = csFeq(:,[2 end-1]);
clFeq([1 end],:) = clFeq([2 end-1],:);
clFeq(:,[1 end]) = clFeq(:,[2 end-1]);

[fsSiq,csSiq,clSiq] = equilibrium(SOL.T,CHM.cSi,SOL.Pt,CHM.TSi1,CHM.TSi2,CHM.cphsSi1,CHM.cphsSi2,...
                                  CHM.perTSi,CHM.perCsSi,CHM.perClSi,CHM.clap,CHM.PhDgSi);
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
CHM.GFe   = alpha.*CHM.GFe + (1-alpha).*((CHM.XFe.*fsFeq-CHM.FsFe)./(4.*NUM.dt));
advn_FFe  = - advect(CHM.FsFe(inz,inx),UsFe(inz,:),WsFe(:,inx),NUM.h,{ADVN,''},[1,2],BCA);
dFFedt    = advn_FFe + CHM.GFe(inz,inx);

CHM.GSi   = alpha.*CHM.GSi + (1-alpha).*((CHM.XSi.*fsSiq-CHM.FsSi)./(4.*NUM.dt));
advn_FSi  = - advect(CHM.FsSi(inz,inx),UsSi(inz,:),WsSi(:,inx),NUM.h,{ADVN,''},[1,2],BCA);
dFSidt    = advn_FSi + CHM.GSi(inz,inx);                                       % total rate of change

CHM.FsFe(inz,inx) = FFeo(inz,inx) + (NUM.theta.*dFFedt + (1-NUM.theta).*dFFedto).*NUM.dt;
CHM.FsSi(inz,inx) = FSio(inz,inx) + (NUM.theta.*dFSidt + (1-NUM.theta).*dFSidto).*NUM.dt;

CHM.FsFe([1 end],:) = CHM.FsFe([2 end-1],:);  CHM.FsFe(:,[1 end]) = CHM.FsFe(:,[2 end-1]);  % apply boundary conditions
CHM.FsSi([1 end],:) = CHM.FsSi([2 end-1],:);  CHM.FsSi(:,[1 end]) = CHM.FsSi(:,[2 end-1]);  % apply boundary conditions

CHM.FsFe = max(0, CHM.FsFe );
CHM.FsSi = max(0, CHM.FsSi );

% update phase fractions [wt]
CHM.fsFe(hasFe) = max(0,min(1, CHM.FsFe(hasFe)./max(TINY,CHM.XFe(hasFe)) ));
CHM.fsSi(hasSi) = max(0,min(1, CHM.FsSi(hasSi)./max(TINY,CHM.XSi(hasSi)) ));

CHM.flFe(hasFe) = 1-CHM.fsFe(hasFe); 
CHM.flSi(hasSi) = 1-CHM.fsSi(hasSi);


%% update phase compositions
KcFe = csFeq./clFeq;
CHM.clFe(hasFe) = CHM.cFe(hasFe)./(CHM.flFe(hasFe) + CHM.fsFe(hasFe).*KcFe(hasFe));
CHM.csFe(hasFe) = CHM.cFe(hasFe)./(CHM.flFe(hasFe)./KcFe(hasFe) + CHM.fsFe(hasFe));

KcSi = csSiq./clSiq;
CHM.clSi(hasSi) = CHM.cSi(hasSi)./(CHM.flSi(hasSi) + CHM.fsSi(hasSi).*KcSi(hasSi));
CHM.csSi(hasSi) = CHM.cSi(hasSi)./(CHM.flSi(hasSi)./KcSi(hasSi) + CHM.fsSi(hasSi));


%% get residual of thermochemical equations from iterative update
normT   = norm(SOL.T    - Ti   ,2)./(norm(SOL.T   ,2)+TINY);
normxFe = norm(CHM.xFe  - xFei ,2)./(norm(CHM.xFe ,2)+TINY);
normcFe = norm(CHM.cFe  - cFei ,2)./(norm(CHM.cFe ,2)+TINY);
normcSi = norm(CHM.cSi  - cSii ,2)./(norm(CHM.cSi ,2)+TINY);
normfFe = norm(CHM.flFe - flFei,2)./(norm(CHM.flFe,2)+TINY);
normfSi = norm(CHM.flSi - flSii,2)./(norm(CHM.flSi,2)+TINY);
resnorm_TC = normT + normxFe + normcSi + normcFe + normfFe + normfSi;