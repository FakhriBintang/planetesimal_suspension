% planetesimal: solve thermo-chemical equations

if NUM.step>0
    NUM.step;
end

% store previous iteration
Ti    = SOL.T;
fFeli = CHM.fFel;
fSili = CHM.fSil;


%% update entropy and Temperature
advn_S  = -advection(MAT.rho.*CHM.xSi.*CHM.fSis.*SOL.sSis ,UsSi,WsSi,NUM.h,NUM.h,NUM.ADVN,'flx')...
        -  advection(MAT.rho.*CHM.xFe.*CHM.fFes.*SOL.sFes ,UsFe,WsFe,NUM.h,NUM.h,NUM.ADVN,'flx')...
        -  advection(MAT.rho.*CHM.xSi.*CHM.fSil.*SOL.sSil ,UlSi,WlSi,NUM.h,NUM.h,NUM.ADVN,'flx')...
        -  advection(MAT.rho.*CHM.xFe.*CHM.fFel.*SOL.sFel ,UlFe,WlFe,NUM.h,NUM.h,NUM.ADVN,'flx');
% enthalpy on the outside
qSz    = - (MAT.ks(1:end-1,:)+MAT.ks(2:end,:))./2 .* ddz(SOL.T,NUM.h);                     % heat diffusion z-flux
qSx    = - (MAT.ks(:,1:end-1)+MAT.ks(:,2:end))./2 .* ddx(SOL.T,NUM.h);                     % heat diffusion x-flux
diff_S(2:end-1,2:end-1) = (- ddz(qSz(:,2:end-1),NUM.h)  ...                    % heat diffusion
                           - ddx(qSx(2:end-1,:),NUM.h));


diss_T = zeros(size(SOL.T));
diss_T(2:end-1,2:end-1) = EntProd ./SOL.T(2:end-1,2:end-1);

dSdt   = advn_S + diff_S + diss_T; % + MAT.Hr;

if NUM.step>0
    % update solution
    SOL.S = So + (NUM.theta.*dSdt + (1-NUM.theta).*dSdto) .* NUM.dt;

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
adv_XFe     = advection(MAT.rho.*CHM.xFe.*CHM.fFes,UsFe,WsFe,NUM.h,NUM.h,NUM.ADVN,'flx')...
            + advection(MAT.rho.*CHM.xFe.*CHM.fFel,UlFe,WlFe,NUM.h,NUM.h,NUM.ADVN,'flx');
dXdt        = - adv_XFe;

if NUM.step>0 && any(CHM.xFe(:)>0) && any(CHM.xFe(:)<1)
    % update solution
    CHM.XFe = XFeo + (NUM.theta.*dXdt + (1-NUM.theta).*dXdto) .* NUM.dt;

    % apply boundaries
    CHM.XFe([1 end],:) = CHM.XFe([2 end-1],:);  CHM.XFe(:,[1 end]) = CHM.XFe(:,[2 end-1]);
else
    CHM.XFe = CHM.xFe.*MAT.rho;
end

% update fertile chemical components
% equations 8a, 8b
advn_CSi  = - advection(MAT.rho.*CHM.xSi.*CHM.fSis.*CHM.csSi,UsSi,WsSi,NUM.h,NUM.h,NUM.ADVN,'flx')...
            - advection(MAT.rho.*CHM.xSi.*CHM.fSil.*CHM.clSi,UlSi,WlSi,NUM.h,NUM.h,NUM.ADVN,'flx');

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

dCSidt    = advn_CSi + diff_CSi;

advn_CFe  = - advection(MAT.rho.*CHM.xFe.*CHM.fFes.*CHM.csFe,UsFe,WsFe,NUM.h,NUM.h,NUM.ADVN,'flx')...
            - advection(MAT.rho.*CHM.xFe.*CHM.fFel.*CHM.clFe,UlFe,WlFe,NUM.h,NUM.h,NUM.ADVN,'flx');

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

dCFedt      = advn_CFe + diff_CFe;

if NUM.step>0
    % update solution
    CHM.CSi = CSio + (NUM.theta.*dCSidt + (1-NUM.theta).*dCSidto) .* NUM.dt;
    CHM.CFe = CFeo + (NUM.theta.*dCFedt + (1-NUM.theta).*dCFedto) .* NUM.dt;

    % apply boundaries
    CHM.CSi([1 end],:) = CHM.CSi([2 end-1],:);  CHM.CSi(:,[1 end]) = CHM.CSi(:,[2 end-1]);
    CHM.CFe([1 end],:) = CHM.CFe([2 end-1],:);  CHM.CFe(:,[1 end]) = CHM.CFe(:,[2 end-1]); 
end

if NUM.step>0
    CHM.XSi = MAT.rho - CHM.XFe;
    CHM.xFe = CHM.XFe./MAT.rho;
    CHM.xSi = 1-CHM.xFe;
    CHM.cSi = CHM.CSi./(CHM.xSi+TINY)./MAT.rho ;
    CHM.cFe = CHM.CFe./(CHM.xFe+TINY)./MAT.rho ;
    SOL.T    = SOL.T0.*exp(SOL.S./MAT.rho./PHY.Cp ...
             + PHY.aT.*(SOL.Pt - P0)./rhoRef./PHY.Cp...
             - MAT.Ds./PHY.Cp);
else
    SOL.S   = MAT.rho.*(PHY.Cp.*log(SOL.T/SOL.T0) - PHY.aT./rhoRef.*(SOL.Pt-P0) + MAT.Ds);
    SOL.H   = SOL.T.*(MAT.rho.*(PHY.Cp + MAT.Ds));
    CHM.XFe = MAT.rho.*CHM.xFe; CHM.XSi = MAT.rho.*CHM.xSi;
    CHM.XSi = MAT.rho - CHM.XFe;
    CHM.CFe = MAT.rho.*CHM.cFe.*CHM.xFe;
    CHM.CSi = MAT.rho.*CHM.cSi.*CHM.xSi;
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
    CHM.GFe     = alpha.*CHM.GFe + (1-alpha).*((fFelq-CHM.fFel).*CHM.XFe./max(4.*NUM.dt,CHM.tau_r));

    advn_FFe    = - advection(MAT.rho.*CHM.xFe.*CHM.fFel,UlFe,WlFe,NUM.h,NUM.h,NUM.ADVN,'flx');            % get advection term

    dFFedt      = advn_FFe + CHM.GFe;

    CHM.GSi     = alpha.*CHM.GSi + (1-alpha).*((fSilq-CHM.fSil).*CHM.XSi./max(4.*NUM.dt,CHM.tau_r));

    advn_FSi    = - advection(MAT.rho.*CHM.xSi.*CHM.fSil,UlSi,WlSi,NUM.h,NUM.h,NUM.ADVN,'flx');            % get advection term

    dFSidt      = advn_FSi + CHM.GSi;                                       % total rate of change



    if NUM.step>0 % update phase 
        FFe = rhoo.*xFeo.*fFelo + (NUM.theta.*dFFedt + (1-NUM.theta).*dFFedto).*NUM.dt; FFe = min(MAT.rho-TINY,max(TINY,FFe));
        FSi = rhoo.*xSio.*fSilo + (NUM.theta.*dFSidt + (1-NUM.theta).*dFSidto).*NUM.dt; FSi = min(MAT.rho-TINY,max(TINY,FSi));
        CHM.fFel = FFe./(CHM.xFe+TINY)./MAT.rho;
        CHM.fSil = FSi./(CHM.xSi+TINY)./MAT.rho;  % explicit update of crystal fractionCHM.fSil = min(1-TINY,max(TINY,CHM.fSil));                             % enforce [0,1] limit
        CHM.fSil = min(1,max(0,CHM.fSil)); CHM.fFel = min(1,max(0,CHM.fFel));
    end

    CHM.fFel([1 end],:) = CHM.fFel([2 end-1],:);                           % apply boundary conditions
    CHM.fFel(:,[1 end]) = CHM.fFel(:,[2 end-1]);

    CHM.fSil([1 end],:) = CHM.fSil([2 end-1],:);                           % apply boundary conditions
    CHM.fSil(:,[1 end]) = CHM.fSil(:,[2 end-1]);
    
    CHM.fFes = 1-CHM.fFel; CHM.fSis = 1-CHM.fSil;

else
    % lag equilibrium phase fractions
    CHM.fFes = alpha.*CHM.fFes + (1-alpha).*fFesq;
    CHM.fSis = alpha.*CHM.fSis + (1-alpha).*fSisq;

    CHM.fFel = 1-CHM.fFes; CHM.fSil = 1-CHM.fSis;

    CHM.GFe = (MAT.rho.*CHM.fFel-rhoo.*fFelo)./NUM.dt + advection(MAT.rho.*CHM.xFe.*CHM.fFel,UlFe,WlFe,NUM.h,NUM.h,NUM.ADVN,'flx');  % reconstruct iron melting rate
    CHM.GSi = (MAT.rho.*CHM.fSil-rhoo.*fSilo)./NUM.dt + advection(MAT.rho.*CHM.xSi.*CHM.fSil,UlSi,WlSi,NUM.h,NUM.h,NUM.ADVN,'flx');  % reconstruct silicate melting rate
end

% update phase compositions
KcFe = csFeq./clFeq;
CHM.clFe = CHM.cFe./(CHM.fFel + CHM.fFes.*KcFe);
CHM.csFe = CHM.cFe./(CHM.fFel./KcFe + CHM.fFes);

KcSi = csSiq./clSiq;
CHM.clSi = CHM.cSi./(CHM.fSil + CHM.fSis.*KcSi);
CHM.csSi = CHM.cSi./(CHM.fSil./KcSi + CHM.fSis);

% update entropies
SOL.sFes = SOL.S./MAT.rho -MAT.Ds;
SOL.sSis = SOL.S./MAT.rho -MAT.Ds;
SOL.sFel = SOL.sFes + CHM.dEntrFe;
SOL.sSil = SOL.sSis + CHM.dEntrSi;
sumS = MAT.rho.*(CHM.xSi.*CHM.fSis.*SOL.sSis +CHM.xSi.*CHM.fSil.*SOL.sSil + CHM.xFe.*CHM.fFes.*SOL.sFes +CHM.xFe.*CHM.fFel.*SOL.sFel);

% get residual of thermochemical equations from iterative update
resnorm_TC  = norm(SOL.T    - Ti   ,2)./norm(SOL.T   ,2) ...
            + norm(CHM.fFel - fFeli,2)./norm(CHM.fFel,2) ...
            + norm(CHM.fSil - fSili,2)./norm(CHM.fSil,2);