% planetesimal: solve thermo-chemical equations

if step>0
    step;
end

% store previous iteration
Ti    = T;
xFei  = xFe;
xSii  = xSi;
cFei  = cFe;
cSii  = cSi;
flFei = flFe;
flSii = flSi;


%% update phase entropies
slFe  = (S - FsFe.*dEntrFe - FsSi.*dEntrSi)./rho;
slSi  = slFe;
ssFe  = slFe + dEntrFe;
ssSi  = slSi + dEntrSi;


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

dSdt   = advn_S + diff_S + Hr + diss_T;

% update solution
S(inz,inx) = So(inz,inx) + (theta.*dSdt + (1-theta).*dSdto) .* dt;

% apply boundaries
Ds = xFe.*fsFe.*dEntrFe + xSi.*fsSi.*dEntrSi;
switch BCTTop
    case 'isothermal'
        S(1,:) = rho(1,:).*(Cp.*log(T0./T0)+aT./rhoRef.*(Pt(1,:)-P0) + Ds(1,:));
    case 'insulating'
        S(1,:) = S(2,:);
    case 'flux'
end
switch BCTBot
    case 'isothermal'
        S(end,:) = rho(end,:).*(Cp.*log(T1./T0)+aT./rhoRef.*(Pt(end,:)-P0) + Ds(end,:));
    case 'insulating'
        S(end,:) = S(end-1,:);
end
switch BCTSides
    case 'isothermal'
        S(:,[1 end]) = rho(:,[1 end]).*(Cp.*log(T0./T0)+aT./rhoRef.*(Pt(:,[1 end])-P0) + Ds(:,[1 end]));
    case 'insulating'
        S(:,[1 end]) = S(:,[2 end-1]);
end

% update temperature
T   = T0.*exp((S - FsFe.*dEntrFe - FsSi.*dEntrSi)./rho./Cp ...
        + aT.*(Pt - P0)./rhoRef./Cp);


%% update system fractions, only one system needs to be solved as SUM_i(X_i) = 1
% equation 7

% % update system indices
% hasFe   = xFe         >0;
% hasSi   = xSi         >0;
% hasFein = xFe(inz,inx)>0;
% hasSiin = xSi(inz,inx)>0;
% 
% % if any(xFe(:)>0 & xFe(:)<1)
%     dXFedt = - advect(FsFe(inz,inx),UsFe(inz,:),WsFe(:,inx),h,{ADVN,''},[1,2],BCA) ...
%              - advect(FlFe(inz,inx),UlFe(inz,:),WlFe(:,inx),h,{ADVN,''},[1,2],BCA);
% 
%     % update solution
%     XFe(inz,inx) = XFeo(inz,inx) + (theta.*dXFedt + (1-theta).*dXFedto) .* dt;
% 
%     dXSidt = - advect(FsSi(inz,inx),UsSi(inz,:),WsSi(:,inx),h,{ADVN,''},[1,2],BCA) ...
%              - advect(FlSi(inz,inx),UlSi(inz,:),WlSi(:,inx),h,{ADVN,''},[1,2],BCA);
%     
%     % update solution
%     XSi(inz,inx) = XSio(inz,inx) + (theta.*dXSidt + (1-theta).*dXSidto) .* dt;
% 
%     XFe(~hasFe) = rho(~hasFe) - XSi(~hasFe);
%     XSi( hasFe) = rho( hasFe) - XFe( hasFe);
% 
%     % apply boundaries
%     XFe([1 end],:) = XFe([2 end-1],:);  XFe(:,[1 end]) = XFe(:,[2 end-1]);
%     XSi([1 end],:) = XSi([2 end-1],:);  XSi(:,[1 end]) = XSi(:,[2 end-1]);
% 
%     % enforce 0,rho limits
%     XFe = max(0, XFe ); 
%     XSi = max(0, XSi ); 
% 
% % else
% %     XFe = xFe.*rho;
% %     XSi = rho - XFe;
% % end
% 
% % update system fractions
% xFe = max(0,min(1, XFe./rho ));
% xSi = max(0,min(1, XSi./rho ));
% 
% hasFe   = xFe         >0;
% hasSi   = xSi         >0;
% hasFein = xFe(inz,inx)>0;
% hasSiin = xSi(inz,inx)>0;


%% alternate method for X conservation
if any(xFe(:)>0 & xFe(:)<1)
    dXFedt = - advect(FsFe(inz,inx),UsFe(inz,:),WsFe(:,inx),h,{ADVN,''},[1,2],BCA) ...
             - advect(FlFe(inz,inx),UlFe(inz,:),WlFe(:,inx),h,{ADVN,''},[1,2],BCA);

    % update solution
    XFe(inz,inx) = XFeo(inz,inx) + (theta.*dXFedt + (1-theta).*dXFedto) .* dt;

    dXSidt = - advect(FsSi(inz,inx),UsSi(inz,:),WsSi(:,inx),h,{ADVN,''},[1,2],BCA) ...
             - advect(FlSi(inz,inx),UlSi(inz,:),WlSi(:,inx),h,{ADVN,''},[1,2],BCA);
    
    % update solution
    XSi(inz,inx) = XSio(inz,inx) + (theta.*dXSidt + (1-theta).*dXSidto) .* dt;

%     XSi = rho - XFe;

    % apply boundaries
    XFe([1 end],:) = XFe([2 end-1],:);  XFe(:,[1 end]) = XFe(:,[2 end-1]);
    XSi([1 end],:) = XSi([2 end-1],:);  XSi(:,[1 end]) = XSi(:,[2 end-1]);

    % enforce 0,rho limits
    XFe = max(0, XFe ); 
    XSi = max(0, XSi ); 

else
    XFe = xFe.*rho;
    XSi = rho - XFe;
end

% update system fractions
xFe = max(0,min(1, XFe./rho ));
xSi = max(0,min(1, XSi./rho ));

hasFe   = xFe>0 & xSi<1;
hasSi   = xSi>0 & xFe<1;

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
CSi(inz,inx) = CSio(inz,inx) + (theta.*dCSidt + (1-theta).*dCSidto) .* dt;
CFe(inz,inx) = CFeo(inz,inx) + (theta.*dCFedt + (1-theta).*dCFedto) .* dt;


% cheat a little bit again and force C_i = X_i*c_i(t-1) when below the
% solidus or above the liquidus
CSi(~hassolSi) = XSi(~hassolSi) .*cSio(~hassolSi);
CFe(~hassolFe) = XFe(~hassolFe) .*cFeo(~hassolFe);
CSi(~hasliqSi) = XSi(~hasliqSi) .*cSio(~hasliqSi);
CFe(~hasliqFe) = XFe(~hasliqFe) .*cFeo(~hasliqFe);

% apply boundaries
CSi([1 end],:) = CSi([2 end-1],:);  CSi(:,[1 end]) = CSi(:,[2 end-1]);
CFe([1 end],:) = CFe([2 end-1],:);  CFe(:,[1 end]) = CFe(:,[2 end-1]);

% % enforce 0,rho limits
CFe = max(0,CFe);
CSi = max(0,CSi);

% update chemical composition
cSi(hasSi) = CSi(hasSi)./XSi(hasSi);
cFe(hasFe) = CFe(hasFe)./XFe(hasFe);


%% update local phase equilibrium
[fsFeq,csFeq,clFeq] = equilibrium(T,cFe,Pt,TFe1,TFe2,cphsFe1,cphsFe2,...
                                  perTFe,perCsFe,perClFe,clap,PhDgFe);
% apply boundaries
fsFeq([1 end],:) = fsFeq([2 end-1],:);
fsFeq(:,[1 end]) = fsFeq(:,[2 end-1]);
csFeq([1 end],:) = csFeq([2 end-1],:);
csFeq(:,[1 end]) = csFeq(:,[2 end-1]);
clFeq([1 end],:) = clFeq([2 end-1],:);
clFeq(:,[1 end]) = clFeq(:,[2 end-1]);

[fsSiq,csSiq,clSiq] = equilibrium(T,cSi,Pt,TSi1,TSi2,cphsSi1,cphsSi2,...
                                  perTSi,perCsSi,perClSi,clap,PhDgSi);
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
GFe   = alpha.*GFe + (1-alpha).*((XFe.*fsFeq-FsFe)./(4.*dt));
advn_FFe  = - advect(FsFe(inz,inx),UsFe(inz,:),WsFe(:,inx),h,{ADVN,''},[1,2],BCA);
dFFedt    = advn_FFe + GFe(inz,inx);

GSi   = alpha.*GSi + (1-alpha).*((XSi.*fsSiq-FsSi)./(4.*dt));
advn_FSi  = - advect(FsSi(inz,inx),UsSi(inz,:),WsSi(:,inx),h,{ADVN,''},[1,2],BCA);
dFSidt    = advn_FSi + GSi(inz,inx);                                       % total rate of change

FsFe(inz,inx) = FFeo(inz,inx) + (theta.*dFFedt + (1-theta).*dFFedto).*dt;
FsSi(inz,inx) = FSio(inz,inx) + (theta.*dFSidt + (1-theta).*dFSidto).*dt;

FsFe([1 end],:) = FsFe([2 end-1],:);  FsFe(:,[1 end]) = FsFe(:,[2 end-1]);  % apply boundary conditions
FsSi([1 end],:) = FsSi([2 end-1],:);  FsSi(:,[1 end]) = FsSi(:,[2 end-1]);  % apply boundary conditions

FsFe = max(0, FsFe );
FsSi = max(0, FsSi );

FlFe = XFe - FsFe;
FlSi = XSi - FsSi;

% update phase fractions [wt]
fsFe(hasFe) = max(0,min(1, FsFe(hasFe)./max(TINY,XFe(hasFe)) ));
fsSi(hasSi) = max(0,min(1, FsSi(hasSi)./max(TINY,XSi(hasSi)) ));

flFe(hasFe) = 1-fsFe(hasFe); 
flSi(hasSi) = 1-fsSi(hasSi);

%detect where is fully molten
hassolSi = flSi<1;
hassolFe = flFe<1;

hasliqSi = fsSi<1;
hasliqFe = fsFe<1;

%% update phase compositions
KcFe = csFeq./clFeq;
clFe(hasFe) = cFe(hasFe)./(flFe(hasFe) + fsFe(hasFe).*KcFe(hasFe));
csFe(hasFe) = cFe(hasFe)./(flFe(hasFe)./KcFe(hasFe) + fsFe(hasFe));

KcSi = csSiq./clSiq;
clSi(hasSi) = cSi(hasSi)./(flSi(hasSi) + fsSi(hasSi).*KcSi(hasSi));
csSi(hasSi) = cSi(hasSi)./(flSi(hasSi)./KcSi(hasSi) + fsSi(hasSi));


%% get residual of thermochemical equations from iterative update
normT   = norm(T    - Ti   ,2)./(norm(T   ,2)+TINY);
normxFe = norm(xFe  - xFei ,2)./(norm(xFe ,2)+TINY);
normcFe = norm(cFe  - cFei ,2)./(norm(cFe ,2)+TINY);
normcSi = norm(cSi  - cSii ,2)./(norm(cSi ,2)+TINY);
normfFe = norm(flFe - flFei,2)./(norm(flFe,2)+TINY);
normfSi = norm(flSi - flSii,2)./(norm(flSi,2)+TINY);
resnorm_TC = normT + normxFe + normcSi + normcFe + normfFe + normfSi;