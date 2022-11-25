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
fFeli = CHM.fFel;
fSili = CHM.fSil;
rhoi  = MAT.rho;
XSii  = CHM.XSi;
XFei  = CHM.XFe;


%% update entropy and Temperature
advn_S  = - advect(MAT.rho(inz,inx).*CHM.xSi(inz,inx).*CHM.fSis(inz,inx).*SOL.sSis(inz,inx),UsSi(inz,:),WsSi(:,inx),NUM.h,{ADVN,''},[1,2],BCA) ...  % heat advection
          - advect(MAT.rho(inz,inx).*CHM.xFe(inz,inx).*CHM.fFes(inz,inx).*SOL.sFes(inz,inx),UsFe(inz,:),WsFe(:,inx),NUM.h,{ADVN,''},[1,2],BCA) ...
          - advect(MAT.rho(inz,inx).*CHM.xSi(inz,inx).*CHM.fSil(inz,inx).*SOL.sSil(inz,inx),UlSi(inz,:),WlSi(:,inx),NUM.h,{ADVN,''},[1,2],BCA) ...  % heat advection
          - advect(MAT.rho(inz,inx).*CHM.xFe(inz,inx).*CHM.fFel(inz,inx).*SOL.sFel(inz,inx),UlFe(inz,:),WlFe(:,inx),NUM.h,{ADVN,''},[1,2],BCA);

% enthalpy on the outside
qSz    = - (MAT.ks(1:end-1,:)+MAT.ks(2:end,:))./2 .* ddz(SOL.T,NUM.h);                     % heat diffusion z-flux
qSx    = - (MAT.ks(:,1:end-1)+MAT.ks(:,2:end))./2 .* ddx(SOL.T,NUM.h);                     % heat diffusion x-flux
diff_S =(- ddz(qSz(:,2:end-1),NUM.h)  ...                    % heat diffusion
         - ddx(qSx(2:end-1,:),NUM.h));


% diss_T = zeros(size(SOL.T));
diss_T = EntProd ./SOL.T(2:end-1,2:end-1);

dSdt   = advn_S + diff_S + diss_T; % + MAT.Hr;

if NUM.step>0
    % update solution
    SOL.S(inz,inx) = So(inz,inx) + (NUM.theta.*dSdt + (1-NUM.theta).*dSdto) .* NUM.dt;

    %apply boundaries
    switch SOL.BCTTop
        case 'isothermal'
            SOL.S(1,:)      = MAT.rho(1,:)  .*(PHY.Cp.*log(SOL.T0./SOL.T0)+PHY.aT./rhoRef.*(SOL.Pt(1,:)  -P0) + MAT.Ds(1,:));
        case 'insulating'
            SOL.S(1,:)      =  SOL.S(2,:);
        case 'flux'
    end
    switch SOL.BCTBot
        case 'isothermal'
            SOL.S(end,:)    = MAT.rho(end,:).*(PHY.Cp.*log(SOL.T1./SOL.T0)+PHY.aT./rhoRef.*(SOL.Pt(end,:)-P0) + MAT.Ds(end,:));                   
        case 'insulating'
            SOL.S(end,:)    =  SOL.S(end-1,:);
    end
    switch SOL.BCTSides
        case 'isothermal'
            SOL.S(:,[1 end])  =  MAT.rho(:,[1 end]).*(PHY.Cp.*log(SOL.T0./SOL.T0)+PHY.aT./rhoRef.*(SOL.Pt(:,[1 end])-P0) + MAT.Ds(:,[1 end]));
        case 'insulating'
            SOL.S(:,[1 end])  =  SOL.S(:,[2 end-1]);
    end

end


%% update composition

% update system, only one system needs to be solved as SUM_i(X_i) = 1
% equation 7

if any(CHM.xFe(:)>0 & CHM.xFe(:)<1)

adv_XFe     = - advect(MAT.rho(inz,inx).*CHM.xFe(inz,inx).*CHM.fFes(inz,inx),UsFe(inz,:),WsFe(:,inx),NUM.h,{ADVN,''},[1,2],BCA) ...
              - advect(MAT.rho(inz,inx).*CHM.xFe(inz,inx).*CHM.fFel(inz,inx),UlFe(inz,:),WlFe(:,inx),NUM.h,{ADVN,''},[1,2],BCA);
dXFedt        = adv_XFe;
dXSidt      = - advect(MAT.rho(inz,inx).*CHM.xSi(inz,inx).*CHM.fSis(inz,inx),UsSi(inz,:),WsSi(:,inx),NUM.h,{ADVN,''},[1,2],BCA) ...
              - advect(MAT.rho(inz,inx).*CHM.xSi(inz,inx).*CHM.fSil(inz,inx),UlSi(inz,:),WlSi(:,inx),NUM.h,{ADVN,''},[1,2],BCA);
if NUM.step>0 && any(CHM.xFe(:)>0) && any(CHM.xFe(:)<1)
    % update solution
    XFe1(inz,inx) = XFeo(inz,inx) + (NUM.theta.*dXFedt + (1-NUM.theta).*dXFedto) .* NUM.dt;
    % apply boundaries
    XFe1([1 end],:) = XFe1([2 end-1],:);  XFe1(:,[1 end]) = XFe1(:,[2 end-1]);
    % enforce 0,rho limits
    XFe1 = min(MAT.rho,max(0,XFe1));   

        % update solution
    XSi1(inz,inx) = XSio(inz,inx) + (NUM.theta.*dXSidt + (1-NUM.theta).*dXSidto) .* NUM.dt;
    % apply boundaries
    XSi1([1 end],:) = XSi1([2 end-1],:);  XSi1(:,[1 end]) = XSi1(:,[2 end-1]);
    % enforce 0,rho limits
    XSi1 = min(MAT.rho,max(0,XSi1)); 

    % normalise system densities
    XFenorm = XFe1./(XFe1+XSi1); XSinorm = XSi1./(XFe1+XSi1); 

    CHM.XFe = MAT.rho.*XFenorm; CHM.XSi = MAT.rho.*XSinorm;
else
    CHM.XFe = CHM.xFe.*MAT.rho;
    CHM.XSi = MAT.rho - CHM.XFe;
end
end

% update fertile chemical components
% equations 8a, 8b

advn_CSi    =-advect (MAT.rho(inz,inx).*CHM.xSi(inz,inx).*CHM.fSis(inz,inx).*CHM.csSi(inz,inx),UsSi(inz,:),WsSi(:,inx),NUM.h,{ADVN,''},[1,2],BCA) ...
             -advect (MAT.rho(inz,inx).*CHM.xSi(inz,inx).*CHM.fSil(inz,inx).*CHM.clSi(inz,inx),UlSi(inz,:),WlSi(:,inx),NUM.h,{ADVN,''},[1,2],BCA);

qCSiz     = - PHY.kC.* (MAT.rho (1:end-1,:)+MAT.rho (2:end,:))/2 ...
          .* (CHM.xSi (1:end-1,:)+CHM.xSi (2:end,:))/2 ...
          .* (CHM.fSil(1:end-1,:)+CHM.fSil(2:end,:))/2 ...
          .* ddz(CHM.cSi,NUM.h);                                 % chemical diffusion z-flux
qCSix     = - PHY.kC.* (MAT.rho (:,1:end-1)+MAT.rho (:,2:end))/2 ...
          .* (CHM.xSi (:,1:end-1)+CHM.xSi (:,2:end))/2 ...
          .* (CHM.fSil(:,1:end-1)+CHM.fSil(:,2:end))/2 ...
          .* ddx(CHM.cSi,NUM.h);
diff_CSi = (- ddz(qCSiz(:,2:end-1),NUM.h) ...             % chemical diffusion
            - ddx(qCSix(2:end-1,:),NUM.h));


dCSidt    = advn_CSi + diff_CSi;

advn_CFe    =-advect (MAT.rho(inz,inx).*CHM.xFe(inz,inx).*CHM.fFes(inz,inx).*CHM.csFe(inz,inx),UsFe(inz,:),WsFe(:,inx),NUM.h,{ADVN,''},[1,2],BCA) ...
             -advect (MAT.rho(inz,inx).*CHM.xFe(inz,inx).*CHM.fFel(inz,inx).*CHM.clFe(inz,inx),UlFe(inz,:),WlFe(:,inx),NUM.h,{ADVN,''},[1,2],BCA);

qCFez     = - PHY.kC.* (MAT.rho (1:end-1,:)+MAT.rho (2:end,:))/2 ...
          .* (CHM.xFe (1:end-1,:)+CHM.xFe (2:end,:))/2 ...
          .* (CHM.fFel(1:end-1,:)+CHM.fFel(2:end,:))/2 ...
          .* ddz(CHM.cFe,NUM.h);                                 % chemical diffusion z-flux
qCFex     = - PHY.kC.* (MAT.rho (:,1:end-1)+MAT.rho (:,2:end))/2 ...
          .* (CHM.xFe (:,1:end-1)+CHM.xFe (:,2:end))/2 ...
          .* (CHM.fFel(:,1:end-1)+CHM.fFel(:,2:end))/2 ...
          .* ddx(CHM.cFe,NUM.h);                                 % chemical diffusion z-flux
diff_CFe  = (- ddz(qCFez(:,2:end-1),NUM.h) ...             % chemical diffusion
             - ddx(qCFex(2:end-1,:),NUM.h)).*0;

dCFedt      = advn_CFe + diff_CFe;

if NUM.step>0
    % update solution
    CHM.CSi(inz,inx) = CSio(inz,inx) + (NUM.theta.*dCSidt + (1-NUM.theta).*dCSidto) .* NUM.dt;
    CHM.CFe(inz,inx) = CFeo(inz,inx) + (NUM.theta.*dCFedt + (1-NUM.theta).*dCFedto) .* NUM.dt;

    % apply boundaries
    CHM.CSi([1 end],:) = CHM.CSi([2 end-1],:);  CHM.CSi(:,[1 end]) = CHM.CSi(:,[2 end-1]);
    CHM.CFe([1 end],:) = CHM.CFe([2 end-1],:);  CHM.CFe(:,[1 end]) = CHM.CFe(:,[2 end-1]); 

    % enforce 0,rho limits
    CHM.CFe = min(CHM.XFe,max(0,CHM.CFe));   
    CHM.CSi = min(CHM.XSi,max(0,CHM.CSi));   
end   




%% update local phase equilibrium
[fFesq,csFeq,clFeq] = equilibrium_single(SOL.T,CHM.cFe,SOL.Pt,CHM.TFe1,CHM.TFe2,CHM.cphsFe1,CHM.cphsFe2,...
    CHM.perTFe,CHM.perCsFe,CHM.perClFe,CHM.clap,CHM.PhDgFe,TINY);

[fSisq,csSiq,clSiq] = equilibrium(SOL.T,CHM.cSi,SOL.Pt,CHM.TSi1,CHM.TSi2,CHM.cphsSi1,CHM.cphsSi2,...
    CHM.perTSi,CHM.perCsSi,CHM.perClSi,CHM.clap,CHM.PhDgSi,TINY);

fFelq = 1-fFesq;
fSilq = 1-fSisq;

% update phase fractions
if RUN.diseq
    CHM.GFe     = alpha.*CHM.GFe + (1-alpha).*((fFesq-CHM.fFes).*XFei./(4.*NUM.dt));
    advn_FFe    = - advect (MAT.rho(inz,inx).*CHM.xFe(inz,inx).*CHM.fFes(inz,inx),UsFe(inz,:),WsFe(:,inx),NUM.h,{ADVN,''},[1,2],BCA);
    dFFedt      = advn_FFe + CHM.GFe(inz,inx);

    CHM.GSi     = alpha.*CHM.GSi + (1-alpha).*((fSisq-CHM.fSis).*XSii./(4.*NUM.dt));
    advn_FSi    = - advect (MAT.rho(inz,inx).*CHM.xSi(inz,inx).*CHM.fSis(inz,inx),UsSi(inz,:),WsSi(:,inx),NUM.h,{ADVN,''},[1,2],BCA);
    dFSidt      = advn_FSi + CHM.GSi(inz,inx);                                       % total rate of change



    if NUM.step>0 % update phase 
        FFe(inz,inx) = FFeo(inz,inx) + (NUM.theta.*dFFedt + (1-NUM.theta).*dFFedto).*NUM.dt; FFe = min(CHM.XFe,max(0,FFe));
        FSi(inz,inx) = FSio(inz,inx) + (NUM.theta.*dFSidt + (1-NUM.theta).*dFSidto).*NUM.dt; FSi = min(CHM.XSi,max(0,FSi));
        FFe([1 end],:) = FFe([2 end-1],:);  FFe(:,[1 end]) = FFe(:,[2 end-1]); % apply boundary conditions
        FSi([1 end],:) = FSi([2 end-1],:);  FSi(:,[1 end]) = FSi(:,[2 end-1]); % apply boundary conditions
        CHM.fFes = FFe./max(TINY,CHM.XFe);
        CHM.fSis = FSi./max(TINY,CHM.XSi);
    end

    CHM.fFes([1 end],:) = CHM.fFes([2 end-1],:);                           % apply boundary conditions
    CHM.fFes(:,[1 end]) = CHM.fFes(:,[2 end-1]);

    CHM.fSis([1 end],:) = CHM.fSis([2 end-1],:);                           % apply boundary conditions
    CHM.fSis(:,[1 end]) = CHM.fSis(:,[2 end-1]);
    
    CHM.fFel = 1-CHM.fFes; CHM.fSil = 1-CHM.fSis;
    CHM.fFel = max(0,min(CHM.fFel,1)); CHM.fSil = max(0,min(CHM.fSil,1));
    CHM.fFes = max(0,min(CHM.fFes,1)); CHM.fSis = max(0,min(CHM.fSis,1));

else
    % equilibrium
    % lag equilibrium phase fractions
    CHM.fFes = alpha.*CHM.fFes + (1-alpha).*fFesq;
    CHM.fSis = alpha.*CHM.fSis + (1-alpha).*fSisq;

    CHM.fFel = 1-CHM.fFes; CHM.fSil = 1-CHM.fSis;

    CHM.GFe = (MAT.rho.*CHM.xFe.*CHM.fFel-rhoo.*xFeo.*fFelo)./NUM.dt ...
            + advect(MAT.rho(inz,inx).*CHM.xFe(inz,inx).*CHM.fFel(inz,inx),UlFe,WlFe,NUM.h,NUM.h,NUM.ADVN,'flx');  % reconstruct iron melting rate
    CHM.GSi = (MAT.rho.*CHM.xSi.*CHM.fSil-rhoo.*xSio.*fSilo)./NUM.dt ...
            + advect(MAT.rho(inz,inx).*CHM.xSi(inz,inx).*CHM.fSil(inz,inx),UlSi,WlSi,NUM.h,NUM.h,NUM.ADVN,'flx');  % reconstruct silicate melting rate
end

% update phase compositions
if any(CHM.xFe(:)>0 & CHM.xFe(:)<1)
KcFe = csFeq./clFeq;
CHM.clFe = CHM.cFe./(CHM.fFel + CHM.fFes.*KcFe);
CHM.csFe = CHM.cFe./(CHM.fFel./KcFe + CHM.fFes);
else
    CHM.clFe(:,:) = 0; CHM.csFe(:,:) = 0;
end

if any(CHM.xSi(:)>0 | CHM.xSi(:)<1)
KcSi = csSiq./clSiq;
CHM.clSi = CHM.cSi./(CHM.fSil + CHM.fSis.*KcSi);
CHM.csSi = CHM.cSi./(CHM.fSil./KcSi + CHM.fSis);
else
    CHM.clSi(:,:) = 0; CHM.csSi(:,:) = 0;
end

% update entropies
SOL.sFes = SOL.S./MAT.rho -MAT.Ds;
SOL.sSis = SOL.S./MAT.rho -MAT.Ds;
SOL.sFel = SOL.sFes + CHM.dEntrFe;
SOL.sSil = SOL.sSis + CHM.dEntrSi;
sumS = MAT.rho.*(CHM.xSi.*CHM.fSis.*SOL.sSis +CHM.xSi.*CHM.fSil.*SOL.sSil...
+ CHM.xFe.*CHM.fFes.*SOL.sFes +CHM.xFe.*CHM.fFel.*SOL.sFel); % testing purposes to compare with S

%% update thermochemical properties from conserved quantities
% update system fractions
CHM.xFe = max(0,min(1,CHM.XFe./MAT.rho));
CHM.xSi = max(0,min(1,CHM.XSi./MAT.rho));
% CHM.xSi = 1 - CHM.xFe;

% update chemical composition
CHM.cSi(CHM.XSi>0) = CHM.CSi(CHM.XSi>0)./CHM.XSi(CHM.XSi>0);
CHM.cFe(CHM.XFe>0) = CHM.CFe(CHM.XFe>0)./CHM.XFe(CHM.XFe>0);
% update temperature
SOL.T   = SOL.T0.*exp(SOL.S./MAT.rho./PHY.Cp ...
    + PHY.aT.*(SOL.Pt - P0)./rhoRef./PHY.Cp...
    - MAT.Ds./PHY.Cp);

%% get residual of thermochemical equations from iterative update
normT   = norm(SOL.T    - Ti   ,2)./(norm(SOL.T   ,2)+TINY);
normxFe = norm(CHM.xFe  - xFei ,2)./(norm(CHM.xFe ,2)+TINY);
normxSi = norm(CHM.xSi  - xSii ,2)./(norm(CHM.xSi ,2)+TINY);
normcFe = norm(CHM.cFe  - cFei ,2)./(norm(CHM.cFe ,2)+TINY);
normcSi = norm(CHM.cSi  - cSii ,2)./(norm(CHM.cSi ,2)+TINY);
normfFe = norm(CHM.fFel - fFeli,2)./(norm(CHM.fFel,2)+TINY);
normfSi = norm(CHM.fSil - fSili,2)./(norm(CHM.fSil,2)+TINY);
normrho = norm(MAT.rho  - rhoi ,2)./(norm(MAT.rho ,2)+TINY);
resnorm_TC  = normT + normxFe + normxSi + normcSi + normcFe + normrho;

% figure(99); if iter==1; clf; else; hold on; end
%             plot(iter,log10(normT),'r.',iter,log10(normxFe),'g.',iter,log10(normxSi),'b.',iter,log10(normcFe),'c.',iter,log10(normcSi),'m.',iter,log10(normfFe),'k.',iter,log10(normfSi),'y.',iter,log10(normrho),'k*','MarkerSize',15,'LineWidth',1.5); box on; axis tight;  
%             drawnow;