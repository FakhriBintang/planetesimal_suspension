t = step+1;

dsumMdto   = dsumMdt;
dsumSdto   = dsumSdt;
dsumXFedto = dsumXFedt;
dsumCFedto = dsumCFedt;
dsumCSidto = dsumCSidt;

HST.time(t) = time;

% record conserved masses at current timestep
HST.sumM  (t)      = sum(sum(rho(2:end-1,2:end-1))) .*h.*h.*1;
HST.sumS  (t)      = sum(sum(  S(2:end-1,2:end-1)+S0(2:end-1,2:end-1))) .*h.*h.*1;
HST.sumXFe(t)      = sum(sum(XFe(2:end-1,2:end-1))) .*h.*h.*1;
HST.sumCFe(t)      = sum(sum(CFe(2:end-1,2:end-1))) .*h.*h.*1;
HST.sumCSi(t)      = sum(sum(CSi(2:end-1,2:end-1))) .*h.*h.*1;

% record expected rates of volume change imposed by variable boundaries
dsumMdt     = sum(rho(2,2:end-1)    .*W(1,2:end-1)  .*h.*1) ...
            - sum(rho(end-1,2:end-1).*W(end,2:end-1).*h.*1) ...
            + sum(rho(2:end-1,2)    .*U(2:end-1,1)  .*h.*1) ...
            - sum(rho(2:end-1,end-1).*U(2:end-1,end).*h.*1);

dsumSdt     = sum(sum(Hr*h*h*1)) + sum(sum(diss_T*h*h*1))...
            + sum(S(2,2:end-1)      .*W(1,2:end-1)  .*h.*1) ...
            - sum(S(end-1,2:end-1)  .*W(end,2:end-1).*h.*1) ...
            + sum(S(2:end-1,2)      .*U(2:end-1,1)  .*h.*1) ...
            - sum(S(2:end-1,end-1)  .*U(2:end-1,end).*h.*1);

dsumXFedt   = sum(XFe(2,2:end-1)    .*W(1,2:end-1)  .*h.*1) ...
            - sum(XFe(end-1,2:end-1).*W(end,2:end-1).*h.*1) ...
            + sum(XFe(2:end-1,2)    .*U(2:end-1,1)  .*h.*1) ...
            - sum(XFe(2:end-1,end-1).*U(2:end-1,end).*h.*1);

dsumCFedt   = sum(CFe(2,2:end-1)    .*W(1,2:end-1)  .*h.*1) ...
            - sum(CFe(end-1,2:end-1).*W(end,2:end-1).*h.*1) ...
            + sum(CFe(2:end-1,2)    .*U(2:end-1,1)  .*h.*1) ...
            - sum(CFe(2:end-1,end-1).*U(2:end-1,end).*h.*1);

dsumCSidt   = sum(CSi(2,2:end-1)    .*W(1,2:end-1)  .*h.*1) ...
            - sum(CSi(end-1,2:end-1).*W(end,2:end-1).*h.*1) ...
            + sum(CSi(2:end-1,2)    .*U(2:end-1,1)  .*h.*1) ...
            - sum(CSi(2:end-1,end-1).*U(2:end-1,end).*h.*1);

if step>0
    HST.dM  (t) = HST.dM  (t-1) + (theta.*dsumMdt   + (1-theta).*dsumMdto  ) .*dt;
    HST.dS  (t) = HST.dS  (t-1) + (theta.*dsumSdt   + (1-theta).*dsumSdto  ) .*dt;
    HST.dXFe(t) = HST.dXFe(t-1) + (theta.*dsumXFedt + (1-theta).*dsumXFedto) .*dt;
    HST.dCFe(t) = HST.dCFe(t-1) + (theta.*dsumCFedt + (1-theta).*dsumCFedto) .*dt;
    HST.dCSi(t) = HST.dCSi(t-1) + (theta.*dsumCSidt + (1-theta).*dsumCSidto) .*dt;
else
    HST.dM  (1) = 0;
    HST.dS  (1) = 0;
    HST.dXFe(1) = 0;
    HST.dCFe(1) = 0;
    HST.dCSi(1) = 0;
end

% Record conservation error (mass - mass change)/initial mass
HST.EM  (t)     = (HST.sumM  (t) - HST.dM  (t))./HST.sumM  (1) - 1;
HST.ES  (t)     = (HST.sumS  (t) - HST.dS  (t))./HST.sumS  (1) - 1;
HST.EXFe(t)     = (HST.sumXFe(t) - HST.dXFe(t))./HST.sumXFe(1) - 1;
HST.ECSi(t)     = (HST.sumCSi(t) - HST.dCSi(t))./HST.sumCSi(1) - 1;
HST.ECFe(t)     = (HST.sumCFe(t) - HST.dCFe(t))./HST.sumCFe(1) - 1;

% record HSTory of solution variables, material properties, and auxiliary fields
HST.T(t,1) = min(min(T(2:end-1,2:end-1)));
HST.T(t,2) = mean(mean(T(2:end-1,2:end-1)));
HST.T(t,3) = max(max(T(2:end-1,2:end-1)));

HST.xFe(t,1) = min(min(xFe(2:end-1,2:end-1)));
HST.xFe(t,2) = mean(mean(xFe(2:end-1,2:end-1)));
HST.xFe(t,3) = max(max(xFe(2:end-1,2:end-1)));

HST.flFe(t,1) = min(min(flFe(2:end-1,2:end-1)));
HST.flFe(t,2) = mean(mean(flFe(2:end-1,2:end-1)));
HST.flFe(t,3) = max(max(flFe(2:end-1,2:end-1)));

HST.flSi(t,1) = min(min(flSi(2:end-1,2:end-1)));
HST.flSi(t,2) = mean(mean(flSi(2:end-1,2:end-1)));
HST.flSi(t,3) = max(max(flSi(2:end-1,2:end-1)));

HST.cFe(t,1) = min(min(cFe(2:end-1,2:end-1)));
HST.cFe(t,2) = mean(mean(cFe(2:end-1,2:end-1)));
HST.cFe(t,3) = max(max(cFe(2:end-1,2:end-1)));

HST.cSi(t,1) = min(min(cSi(2:end-1,2:end-1)));
HST.cSi(t,2) = mean(mean(cSi(2:end-1,2:end-1)));
HST.cSi(t,3) = max(max(cSi(2:end-1,2:end-1)));

HST.phisSi(t,1) = min(min(phisSi(2:end-1,2:end-1)));
HST.phisSi(t,2) = mean(mean(phisSi(2:end-1,2:end-1)));
HST.phisSi(t,3) = max(max(phisSi(2:end-1,2:end-1)));

HST.philSi(t,1) = min(min(philSi(2:end-1,2:end-1)));
HST.philSi(t,2) = mean(mean(philSi(2:end-1,2:end-1)));
HST.philSi(t,3) = max(max(philSi(2:end-1,2:end-1)));

HST.phisFe(t,1) = min(min(phisFe(2:end-1,2:end-1)));
HST.phisFe(t,2) = mean(mean(phisFe(2:end-1,2:end-1)));
HST.phisFe(t,3) = max(max(phisFe(2:end-1,2:end-1)));

HST.philFe(t,1) = min(min(philFe(2:end-1,2:end-1)));
HST.philFe(t,2) = mean(mean(philFe(2:end-1,2:end-1)));
HST.philFe(t,3) = max(max(philFe(2:end-1,2:end-1)));

HST.rho(t,1) = min(min(rho(2:end-1,2:end-1)));
HST.rho(t,2) = mean(mean(rho(2:end-1,2:end-1)));
HST.rho(t,3) = max(max(rho(2:end-1,2:end-1)));

HST.eta(t,1) = min(min(Eta(2:end-1,2:end-1)));
HST.eta(t,2) = geomean(geomean(Eta(2:end-1,2:end-1)));
HST.eta(t,3) = max(max(Eta(2:end-1,2:end-1)));
