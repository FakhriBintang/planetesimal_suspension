t = NUM.step+1;

% record conserved masses at current timestep
CON.Mass  (t)      = sum(sum(MAT.rho(2:end-1,2:end-1))) .*NUM.h.*NUM.h.*1;
CON.sumH  (t)      = sum(sum(SOL.H  (2:end-1,2:end-1))) .*NUM.h.*NUM.h.*1;
CON.sumXFe(t)      = sum(sum(CHM.XFe(2:end-1,2:end-1))) .*NUM.h.*NUM.h.*1;
CON.sumCFe(t)      = sum(sum(CHM.CFe(2:end-1,2:end-1))) .*NUM.h.*NUM.h.*1;
CON.sumCSi(t)      = sum(sum(CHM.CSi(2:end-1,2:end-1))) .*NUM.h.*NUM.h.*1;

% record expected rates of volume change imposed by variable boundaries
dMassdt     = sum(MAT.rho(2,2:end-1)    .*SOL.W(1,2:end-1)  .*NUM.h.*1) ...
            - sum(MAT.rho(end-1,2:end-1).*SOL.W(end,2:end-1).*NUM.h.*1) ...
            + sum(MAT.rho(2:end-1,2)    .*SOL.U(2:end-1,1)  .*NUM.h.*1) ...
            - sum(MAT.rho(2:end-1,end-1).*SOL.U(2:end-1,end).*NUM.h.*1);

dsumHdt     = sum(SOL.H(2,2:end-1)      .*SOL.W(1,2:end-1)  .*NUM.h.*1) ...
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
    CON.dM  (t) = CON.dM  (t-1) + dMassdt   .*NUM.dt;
    CON.dH  (t) = CON.dH  (t-1) + dsumHdt   .*NUM.dt;
    CON.dXFe(t) = CON.dXFe(t-1) + dsumXFedt .*NUM.dt;
    CON.dCFe(t) = CON.dCFe(t-1) + dsumCFedt .*NUM.dt;
    CON.dCSi(t) = CON.dCSi(t-1) + dsumCSidt .*NUM.dt;
else
    CON.dM  (1) = 0;
    CON.dH  (1) = 0;
    CON.dXFe(1) = 0;
    CON.dCFe(1) = 0;
    CON.dCSi(1) = 0;
end

% Record conservation error (mass - mass change)/initial mass
CON.EM  (t)     = (CON.Mass  (t) - CON.dM  (t))./CON.Mass  (1) - 1;
CON.EH  (t)     = (CON.sumH  (t) - CON.dH  (t))./CON.sumH  (1) - 1;
CON.EXFe(t)     = (CON.sumXFe(t) - CON.dXFe(t))./CON.sumXFe(1) - 1;
CON.ECSi(t)     = (CON.sumCSi(t) - CON.dCSi(t))./CON.sumCSi(1) - 1;
CON.ECFe(t)     = (CON.sumCFe(t) - CON.dCFe(t))./CON.sumCFe(1) - 1;
