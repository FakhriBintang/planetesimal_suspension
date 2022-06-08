% planetesimal: solve thermo-chemical equations

% store previous iteration
Ti    = SOL.T;
fFeli = CHM.fFel;
fSili = CHM.fSil;


%% update enthalpy and Temperature
advn_H  = advection(MAT.rho.*CHM.xSi.*CHM.fSis.*SOL.T.* PHY.Cp               ,UsSi,WsSi,NUM.h,NUM.h,NUM.ADVN,'flx')...
        + advection(MAT.rho.*CHM.xFe.*CHM.fFes.*SOL.T.* PHY.Cp               ,UsFe,WsFe,NUM.h,NUM.h,NUM.ADVN,'flx')...
        + advection(MAT.rho.*CHM.xSi.*CHM.fSil.*SOL.T.*(PHY.Cp + CHM.dEntrSi),UlSi,WlSi,NUM.h,NUM.h,NUM.ADVN,'flx')...
        + advection(MAT.rho.*CHM.xFe.*CHM.fFel.*SOL.T.*(PHY.Cp + CHM.dEntrFe),UlFe,WlFe,NUM.h,NUM.h,NUM.ADVN,'flx');

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
            SOL.H(1,:)      =  2*SOL.T0.*(MAT.rho(1,:).*(MAT.Ds(1,:) +PHY.Cp)) - SOL.H(2,:);
        case 'insulating'
            SOL.H(1,:)      =  SOL.H(2,:);
        case 'flux'

    end
    switch SOL.BCTBot
        case 'isothermal'
            SOL.H(end,:)    =  2*SOL.T1.*(MAT.rho(end,:).*(MAT.Ds(end,:) +PHY.Cp)) - SOL.H(end-1,:);
        case 'insulating'
            SOL.H(end,:)    =  SOL.H(end-1,:);
    end    
    switch SOL.BCTSides
        case 'isothermal'
            SOL.H(:,[1 end])  =  SOL.T0.*(MAT.rho(:,[1 end]).*(MAT.Ds(:,[1 end]) + PHY.Cp));
        case 'insulating'
            SOL.H(:,[1 end])  =  SOL.H(:,[2 end-1]);
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
end

if NUM.step>0
    CHM.XSi = MAT.rho - CHM.XFe;
    CHM.xFe = CHM.XFe./MAT.rho;    
    CHM.xSi = 1-CHM.xFe;
    CHM.cSi = CHM.CSi./(CHM.xSi+TINY)./MAT.rho;
    CHM.cFe = CHM.CFe./(CHM.xFe+TINY)./MAT.rho;
    SOL.T   = SOL.H./(MAT.rho.*(PHY.Cp + MAT.Ds));
else
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
    if NUM.step>0; CHM.fFel = (rhoo.*xFeo.*fFelo + (NUM.theta.*dfFedt + (1-NUM.theta).*dfFedto).*NUM.dt)./(CHM.xFe+TINY)./MAT.rho; 
                   CHM.fSil = (rhoo.*xSio.*fSilo + (NUM.theta.*dfSidt + (1-NUM.theta).*dfSidto).*NUM.dt)./(CHM.xSi+TINY)./MAT.rho;end  % explicit update of crystal fraction
    CHM.GFe = alpha.*CHM.GFe + (1-alpha).*((fFelq-CHM.fFel).*MAT.rho./max(4.*NUM.dt,CHM.tau_r));
    
    advn_fFe = advection(MAT.rho.*CHM.xFe.*CHM.fFel,UlFe,WlFe,NUM.h,NUM.h,NUM.ADVN,'flx');            % get advection term
    
    dfFedt   = - advn_fFe + CHM.GFe;            
               
    CHM.fFel = min(1-TINY,max(TINY,CHM.fFel));                             % enforce [0,1] limit
    CHM.fFel([1 end],:) = CHM.fFel([2 end-1],:);                           % apply boundary conditions
    CHM.fFel(:,[1 end]) = CHM.fFel(:,[2 end-1]);

    CHM.GSi = alpha.*CHM.GSi + (1-alpha).*((fSilq-CHM.fSil).*MAT.rho./max(4.*NUM.dt,CHM.tau_r));
    
    advn_fSi = advection(MAT.rho.*CHM.xSi.*CHM.fSil,UlSi,WlSi,NUM.h,NUM.h,NUM.ADVN,'flx');            % get advection term
    
    dfSidt   = - advn_fSi + CHM.GSi;                                       % total rate of change
    
    CHM.fSil = min(1-TINY,max(TINY,CHM.fSil));                             % enforce [0,1] limit
    CHM.fSil([1 end],:) = CHM.fSil([2 end-1],:);                           % apply boundary conditions
    CHM.fSil(:,[1 end]) = CHM.fSil(:,[2 end-1]);
    
    CHM.fFes = 1-CHM.fFel; CHM.fSis = 1-CHM.fSil;

else
    % lag equilibrium phase fractions
    CHM.fFes = alpha.*CHM.fFes + (1-alpha).*fFesq;
    CHM.fSis = alpha.*CHM.fSis + (1-alpha).*fSisq;
    
    CHM.fFel = 1-CHM.fFes; CHM.fSil = 1-CHM.fSis;
    
    CHM.GFe = (MAT.rho.*CHM.fFel-rhoo.*fFelo)./NUM.dt + advection(MAT.rho.*CHM.fFel,UlFe,WlFe,NUM.h,NUM.h,NUM.ADVN,'flx');  % reconstruct iron melting rate
    CHM.GSi = (MAT.rho.*CHM.fSil-rhoo.*fSilo)./NUM.dt + advection(MAT.rho.*CHM.fSil,UlSi,WlSi,NUM.h,NUM.h,NUM.ADVN,'flx');  % reconstruct silicate melting rate
end

% update phase compositions
KcFe = csFeq./clFeq;
CHM.clFe = CHM.cFe./(CHM.fFel + CHM.fFes.*KcFe);
CHM.csFe = CHM.cFe./(CHM.fFel./KcFe + CHM.fFes);

KcSi = csSiq./clSiq;
CHM.clSi = CHM.cSi./(CHM.fSil + CHM.fSis.*KcSi);
CHM.csSi = CHM.cSi./(CHM.fSil./KcSi + CHM.fSis);

% get residual of thermochemical equations from iterative update
resnorm_TC = norm(SOL.T    - Ti   ,2)./norm(SOL.T   ,2) ...
           + norm(CHM.fFel - fFeli,2)./norm(CHM.fFel,2) ...
           + norm(CHM.fSil - fSili,2)./norm(CHM.fSil,2);