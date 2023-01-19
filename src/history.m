t = NUM.step+1;

dsumMdto   = dsumMdt;
dsumSdto   = dsumSdt;
dsumXFedto = dsumXFedt;
dsumCFedto = dsumCFedt;
dsumCSidto = dsumCSidt;

HST.time(t) = NUM.time;

% record conserved masses at current timestep
HST.sumM  (t)      = sum(sum(MAT.rho(2:end-1,2:end-1))) .*NUM.h.*NUM.h.*1;
HST.sumS  (t)      = sum(sum(  SOL.S(2:end-1,2:end-1)+S0(2:end-1,2:end-1))) .*NUM.h.*NUM.h.*1;
HST.sumXFe(t)      = sum(sum(CHM.XFe(2:end-1,2:end-1))) .*NUM.h.*NUM.h.*1;
HST.sumCFe(t)      = sum(sum(CHM.CFe(2:end-1,2:end-1))) .*NUM.h.*NUM.h.*1;
HST.sumCSi(t)      = sum(sum(CHM.CSi(2:end-1,2:end-1))) .*NUM.h.*NUM.h.*1;

% record expected rates of volume change imposed by variable boundaries
dsumMdt     = sum(MAT.rho(2,2:end-1)    .*SOL.W(1,2:end-1)  .*NUM.h.*1) ...
            - sum(MAT.rho(end-1,2:end-1).*SOL.W(end,2:end-1).*NUM.h.*1) ...
            + sum(MAT.rho(2:end-1,2)    .*SOL.U(2:end-1,1)  .*NUM.h.*1) ...
            - sum(MAT.rho(2:end-1,end-1).*SOL.U(2:end-1,end).*NUM.h.*1);

dsumSdt     = sum(sum(MAT.Hr(2:end-1,2:end-1)*NUM.h*NUM.h*1)) ...
            + sum(SOL.S(2,2:end-1)      .*SOL.W(1,2:end-1)  .*NUM.h.*1) ...
            - sum(SOL.S(end-1,2:end-1)  .*SOL.W(end,2:end-1).*NUM.h.*1) ...
            + sum(SOL.S(2:end-1,2)      .*SOL.U(2:end-1,1)  .*NUM.h.*1) ...
            - sum(SOL.S(2:end-1,end-1)  .*SOL.U(2:end-1,end).*NUM.h.*1);

dsumXFedt   = sum(CHM.XFe(2,2:end-1)    .*SOL.W(1,2:end-1)  .*NUM.h.*1) ...
            - sum(CHM.XFe(end-1,2:end-1).*SOL.W(end,2:end-1).*NUM.h.*1) ...
            + sum(CHM.XFe(2:end-1,2)    .*SOL.U(2:end-1,1)  .*NUM.h.*1) ...
            - sum(CHM.XFe(2:end-1,end-1).*SOL.U(2:end-1,end).*NUM.h.*1);

dsumCFedt   = sum(CHM.CFe(2,2:end-1)    .*SOL.W(1,2:end-1)  .*NUM.h.*1) ...
            - sum(CHM.CFe(end-1,2:end-1).*SOL.W(end,2:end-1).*NUM.h.*1) ...
            + sum(CHM.CFe(2:end-1,2)    .*SOL.U(2:end-1,1)  .*NUM.h.*1) ...
            - sum(CHM.CFe(2:end-1,end-1).*SOL.U(2:end-1,end).*NUM.h.*1);

dsumCSidt   = sum(CHM.CSi(2,2:end-1)    .*SOL.W(1,2:end-1)  .*NUM.h.*1) ...
            - sum(CHM.CSi(end-1,2:end-1).*SOL.W(end,2:end-1).*NUM.h.*1) ...
            + sum(CHM.CSi(2:end-1,2)    .*SOL.U(2:end-1,1)  .*NUM.h.*1) ...
            - sum(CHM.CSi(2:end-1,end-1).*SOL.U(2:end-1,end).*NUM.h.*1);

if NUM.step>0
    HST.dM  (t) = HST.dM  (t-1) + (NUM.theta.*dsumMdt   + (1-NUM.theta).*dsumMdto  ) .*NUM.dt;
    HST.dS  (t) = HST.dS  (t-1) + (NUM.theta.*dsumSdt   + (1-NUM.theta).*dsumSdto  ) .*NUM.dt;
    HST.dXFe(t) = HST.dXFe(t-1) + (NUM.theta.*dsumXFedt + (1-NUM.theta).*dsumXFedto) .*NUM.dt;
    HST.dCFe(t) = HST.dCFe(t-1) + (NUM.theta.*dsumCFedt + (1-NUM.theta).*dsumCFedto) .*NUM.dt;
    HST.dCSi(t) = HST.dCSi(t-1) + (NUM.theta.*dsumCSidt + (1-NUM.theta).*dsumCSidto) .*NUM.dt;
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
HST.T(t,1) = min(min(SOL.T(2:end-1,2:end-1)));
HST.T(t,2) = mean(mean(SOL.T(2:end-1,2:end-1)));
HST.T(t,3) = max(max(SOL.T(2:end-1,2:end-1)));

HST.xFe(t,1) = min(min(CHM.xFe(2:end-1,2:end-1)));
HST.xFe(t,2) = mean(mean(CHM.xFe(2:end-1,2:end-1)));
HST.xFe(t,3) = max(max(CHM.xFe(2:end-1,2:end-1)));

HST.flFe(t,1) = min(min(CHM.flFe(2:end-1,2:end-1)));
HST.flFe(t,2) = mean(mean(CHM.flFe(2:end-1,2:end-1)));
HST.flFe(t,3) = max(max(CHM.flFe(2:end-1,2:end-1)));

HST.flSi(t,1) = min(min(CHM.flSi(2:end-1,2:end-1)));
HST.flSi(t,2) = mean(mean(CHM.flSi(2:end-1,2:end-1)));
HST.flSi(t,3) = max(max(CHM.flSi(2:end-1,2:end-1)));

HST.cFe(t,1) = min(min(CHM.cFe(2:end-1,2:end-1)));
HST.cFe(t,2) = mean(mean(CHM.cFe(2:end-1,2:end-1)));
HST.cFe(t,3) = max(max(CHM.cFe(2:end-1,2:end-1)));

HST.cSi(t,1) = min(min(CHM.cSi(2:end-1,2:end-1)));
HST.cSi(t,2) = mean(mean(CHM.cSi(2:end-1,2:end-1)));
HST.cSi(t,3) = max(max(CHM.cSi(2:end-1,2:end-1)));

HST.phisSi(t,1) = min(min(MAT.phisSi(2:end-1,2:end-1)));
HST.phisSi(t,2) = mean(mean(MAT.phisSi(2:end-1,2:end-1)));
HST.phisSi(t,3) = max(max(MAT.phisSi(2:end-1,2:end-1)));

HST.philSi(t,1) = min(min(MAT.philSi(2:end-1,2:end-1)));
HST.philSi(t,2) = mean(mean(MAT.philSi(2:end-1,2:end-1)));
HST.philSi(t,3) = max(max(MAT.philSi(2:end-1,2:end-1)));

HST.phisFe(t,1) = min(min(MAT.phisFe(2:end-1,2:end-1)));
HST.phisFe(t,2) = mean(mean(MAT.phisFe(2:end-1,2:end-1)));
HST.phisFe(t,3) = max(max(MAT.phisFe(2:end-1,2:end-1)));

HST.philFe(t,1) = min(min(MAT.philFe(2:end-1,2:end-1)));
HST.philFe(t,2) = mean(mean(MAT.philFe(2:end-1,2:end-1)));
HST.philFe(t,3) = max(max(MAT.philFe(2:end-1,2:end-1)));

HST.rho(t,1) = min(min(MAT.rho(2:end-1,2:end-1)));
HST.rho(t,2) = mean(mean(MAT.rho(2:end-1,2:end-1)));
HST.rho(t,3) = max(max(MAT.rho(2:end-1,2:end-1)));

HST.eta(t,1) = min(min(MAT.Eta(2:end-1,2:end-1)));
HST.eta(t,2) = geomean(geomean(MAT.Eta(2:end-1,2:end-1)));
HST.eta(t,3) = max(max(MAT.Eta(2:end-1,2:end-1)));
