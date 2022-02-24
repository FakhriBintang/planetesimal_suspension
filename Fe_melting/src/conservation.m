t = NUM.step+1;
% record conserved masses at current timestep
if NUM.step >0
    CON.Mass(t)        = sum(sum(MAT.rhot(2:end-1,2:end-1)))   .*NUM.h.*NUM.h.*1;
    CON.sumH(t)        = sum(sum(SOL.H  (2:end-1,2:end-1)))    .*NUM.h.*NUM.h.*1;
    CON.sumX(t)        = sum(sum(CHM.XFe(2:end-1,2:end-1)))    .*NUM.h.*NUM.h.*1;
    CON.sumCFe(t)      = sum(sum(CHM.CFe(2:end-1,2:end-1)))    .*NUM.h.*NUM.h.*1;
    CON.sumCSi(t)      = sum(sum(CHM.CSi(2:end-1,2:end-1)))    .*NUM.h.*NUM.h.*1;

else 
    CON.Mass(1)        = Mass0;
    CON.sumH(1)        = sumH0;
    CON.sumX(1)        = sumX0;
    CON.sumCFe(1)      = sumCFe0;
    CON.sumCSi(1)      = sumCSi0;
end

% record expected rates of volume change imposed by variable boundaries
dMassdt     = sum(MAT.rhot(2,2:end-1)    .*SOL.W(1,2:end-1)  .*NUM.h.*1) ...
            - sum(MAT.rhot(end-1,2:end-1).*SOL.W(end,2:end-1).*NUM.h.*1) ...
            + sum(MAT.rhot(2:end-1,2)    .*SOL.U(2:end-1,1)  .*NUM.h.*1) ...
            - sum(MAT.rhot(2:end-1,end-1).*SOL.U(2:end-1,end).*NUM.h.*1);

dsumHdt     = sum(SOL.H(2,2:end-1)       .*SOL.W(1,2:end-1)  .*NUM.h.*1) ...
            - sum(SOL.H(end-1,2:end-1)   .*SOL.W(end,2:end-1).*NUM.h.*1) ...
            + sum(SOL.H(2:end-1,2)       .*SOL.U(2:end-1,1)  .*NUM.h.*1) ...
            - sum(SOL.H(2:end-1,end-1)   .*SOL.U(2:end-1,end).*NUM.h.*1);

dsumXdt     = sum(CHM.XFe(2,2:end-1)     .*SOL.W(1,2:end-1)  .*NUM.h.*1) ...
            - sum(CHM.XFe(end-1,2:end-1) .*SOL.W(end,2:end-1).*NUM.h.*1) ...
            + sum(CHM.XFe(2:end-1,2)     .*SOL.U(2:end-1,1)  .*NUM.h.*1) ...
            - sum(CHM.XFe(2:end-1,end-1) .*SOL.U(2:end-1,end).*NUM.h.*1);

dsumCFedt   = sum(CHM.CFe(2,2:end-1)     .*SOL.W(1,2:end-1)  .*NUM.h.*1) ...
            - sum(CHM.CFe(end-1,2:end-1) .*SOL.W(end,2:end-1).*NUM.h.*1) ...
            + sum(CHM.CFe(2:end-1,2)     .*SOL.U(2:end-1,1)  .*NUM.h.*1) ...
            - sum(CHM.CFe(2:end-1,end-1) .*SOL.U(2:end-1,end).*NUM.h.*1);

dsumCSidt   = sum(CHM.CSi(2,2:end-1)     .*SOL.W(1,2:end-1)  .*NUM.h.*1) ...
            - sum(CHM.CSi(end-1,2:end-1) .*SOL.W(end,2:end-1).*NUM.h.*1) ...
            + sum(CHM.CSi(2:end-1,2)     .*SOL.U(2:end-1,1)  .*NUM.h.*1) ...
            - sum(CHM.CSi(2:end-1,end-1) .*SOL.U(2:end-1,end).*NUM.h.*1);

if NUM.step>0
    CON.dM(t)   = CON.dM(t-1)   + dMassdt   .*NUM.dt;
    CON.dH(t)   = CON.dH(t-1)   + dsumHdt   .*NUM.dt;
    CON.dX(t)   = CON.dX(t-1)   + dsumXdt   .*NUM.dt;
    CON.dCFe(t) = CON.dCFe(t-1) + dsumCFedt .*NUM.dt;
    CON.dCSi(t) = CON.dCSi(t-1) + dsumCSidt .*NUM.dt;
else
    CON.dM(1)   = 0;
    CON.dH(1)   = 0;
    CON.dX(1)   = 0;
    CON.dCFe(1) = 0;
    CON.dCSi(1) = 0;
end

% Record conservation error (Mass - mass change)/initial mass, error = 1
% ideal
CON.EM(t)       = (CON.Mass(t)   - CON.dM(t))  ./CON.Mass(1);
CON.EH(t)       = (CON.sumH(t)   - CON.dH(t))  ./CON.sumH(1);
CON.EX(t)       = (CON.sumX(t)   - CON.dX(t))  ./CON.sumX(1);
CON.ECSi(t)     = (CON.sumCSi(t) - CON.dCSi(t))./CON.sumCSi(1);
CON.ECFe(t)     = (CON.sumCFe(t) - CON.dCFe(t))./CON.sumCFe(1);
