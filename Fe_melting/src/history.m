t = NUM.step+1;

HST.time(t) = NUM.time;

% record conserved masses at current timestep
HST.Mass  (t)      = sum(sum(MAT.rho(2:end-1,2:end-1))) .*NUM.h.*NUM.h.*1;
HST.sumH  (t)      = sum(sum(SOL.H  (2:end-1,2:end-1))) .*NUM.h.*NUM.h.*1;
HST.sumXFe(t)      = sum(sum(CHM.XFe(2:end-1,2:end-1))) .*NUM.h.*NUM.h.*1;
HST.sumCFe(t)      = sum(sum(CHM.CFe(2:end-1,2:end-1))) .*NUM.h.*NUM.h.*1;
HST.sumCSi(t)      = sum(sum(CHM.CSi(2:end-1,2:end-1))) .*NUM.h.*NUM.h.*1;

% record expected rates of volume change imposed by variable boundaries
dMassdt     = sum(MAT.rho(2,2:end-1)    .*SOL.W(1,2:end-1)  .*NUM.h.*1) ...
            - sum(MAT.rho(end-1,2:end-1).*SOL.W(end,2:end-1).*NUM.h.*1) ...
            + sum(MAT.rho(2:end-1,2)    .*SOL.U(2:end-1,1)  .*NUM.h.*1) ...
            - sum(MAT.rho(2:end-1,end-1).*SOL.U(2:end-1,end).*NUM.h.*1);

dsumHdt     = sum(sum(MAT.Hr(2:end-1,2:end-1)*NUM.h*NUM.h*1)) ...
            + sum(SOL.H(2,2:end-1)      .*SOL.W(1,2:end-1)  .*NUM.h.*1) ...
            - sum(SOL.H(end-1,2:end-1)  .*SOL.W(end,2:end-1).*NUM.h.*1) ...
            + sum(SOL.H(2:end-1,2)      .*SOL.U(2:end-1,1)  .*NUM.h.*1) ...
            - sum(SOL.H(2:end-1,end-1)  .*SOL.U(2:end-1,end).*NUM.h.*1);

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
    HST.dM  (t) = HST.dM  (t-1) + dMassdt   .*NUM.dt;
    HST.dH  (t) = HST.dH  (t-1) + dsumHdt   .*NUM.dt;
    HST.dXFe(t) = HST.dXFe(t-1) + dsumXFedt .*NUM.dt;
    HST.dCFe(t) = HST.dCFe(t-1) + dsumCFedt .*NUM.dt;
    HST.dCSi(t) = HST.dCSi(t-1) + dsumCSidt .*NUM.dt;
else
    HST.dM  (1) = 0;
    HST.dH  (1) = 0;
    HST.dXFe(1) = 0;
    HST.dCFe(1) = 0;
    HST.dCSi(1) = 0;
end

% Record conservation error (mass - mass change)/initial mass
HST.EM  (t)     = (HST.Mass  (t) - HST.dM  (t))./HST.Mass  (1) - 1;
HST.EH  (t)     = (HST.sumH  (t) - HST.dH  (t))./HST.sumH  (1) - 1;
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

HST.fFel(t,1) = min(min(CHM.fFel(2:end-1,2:end-1)));
HST.fFel(t,2) = mean(mean(CHM.fFel(2:end-1,2:end-1)));
HST.fFel(t,3) = max(max(CHM.fFel(2:end-1,2:end-1)));

HST.fSil(t,1) = min(min(CHM.fSil(2:end-1,2:end-1)));
HST.fSil(t,2) = mean(mean(CHM.fSil(2:end-1,2:end-1)));
HST.fSil(t,3) = max(max(CHM.fSil(2:end-1,2:end-1)));

HST.cFe(t,1) = min(min(CHM.cFe(2:end-1,2:end-1)));
HST.cFe(t,2) = mean(mean(CHM.cFe(2:end-1,2:end-1)));
HST.cFe(t,3) = max(max(CHM.cFe(2:end-1,2:end-1)));

HST.cSi(t,1) = min(min(CHM.cSi(2:end-1,2:end-1)));
HST.cSi(t,2) = mean(mean(CHM.cSi(2:end-1,2:end-1)));
HST.cSi(t,3) = max(max(CHM.cSi(2:end-1,2:end-1)));

HST.phiSis(t,1) = min(min(MAT.phiSis(2:end-1,2:end-1)));
HST.phiSis(t,2) = mean(mean(MAT.phiSis(2:end-1,2:end-1)));
HST.phiSis(t,3) = max(max(MAT.phiSis(2:end-1,2:end-1)));

HST.phiSil(t,1) = min(min(MAT.phiSil(2:end-1,2:end-1)));
HST.phiSil(t,2) = mean(mean(MAT.phiSil(2:end-1,2:end-1)));
HST.phiSil(t,3) = max(max(MAT.phiSil(2:end-1,2:end-1)));

HST.phiFes(t,1) = min(min(MAT.phiFes(2:end-1,2:end-1)));
HST.phiFes(t,2) = mean(mean(MAT.phiFes(2:end-1,2:end-1)));
HST.phiFes(t,3) = max(max(MAT.phiFes(2:end-1,2:end-1)));

HST.phiFel(t,1) = min(min(MAT.phiFel(2:end-1,2:end-1)));
HST.phiFel(t,2) = mean(mean(MAT.phiFel(2:end-1,2:end-1)));
HST.phiFel(t,3) = max(max(MAT.phiFel(2:end-1,2:end-1)));

HST.rho(t,1) = min(min(MAT.rho(2:end-1,2:end-1)));
HST.rho(t,2) = mean(mean(MAT.rho(2:end-1,2:end-1)));
HST.rho(t,3) = max(max(MAT.rho(2:end-1,2:end-1)));

HST.eta(t,1) = min(min(MAT.Eta(2:end-1,2:end-1)));
HST.eta(t,2) = geomean(geomean(MAT.Eta(2:end-1,2:end-1)));
HST.eta(t,3) = max(max(MAT.Eta(2:end-1,2:end-1)));
