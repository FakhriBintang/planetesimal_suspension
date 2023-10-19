t = step+1;

dsumMdtoo  = dsumMdto; dsumMdto = dsumMdt;
dsumSdtoo = dsumSdto; dsumSdto   = dsumSdt;
dsumXFedtoo = dsumXFedto; dsumXFedto = dsumXFedt;
dsumXSidtoo = dsumXSidto; dsumXSidto = dsumXSidt;
dsumCFedtoo = dsumCFedto; dsumCFedto = dsumCFedt;
dsumCSidtoo = dsumCSidto; dsumCSidto = dsumCSidt;
stp = max(1,step+1);
HST.time(t) = time;

% record conserved masses at current timestep
HST.n26Al (t)      = n26Al;
HST.sumM  (t)      = sum(sum(rho(2:end-1,2:end-1))) .*h.*h.*1;
HST.sumS  (t)      = sum(sum(  S(2:end-1,2:end-1)+S0(2:end-1,2:end-1))) .*h.*h.*1;
HST.sumXFe(t)      = sum(sum(XFe(2:end-1,2:end-1))) .*h.*h.*1;
HST.sumXSi(t)      = sum(sum(XSi(2:end-1,2:end-1))) .*h.*h.*1;
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

dsumXSidt   = sum(XSi(2,2:end-1)    .*W(1,2:end-1)  .*h.*1) ...
            - sum(XSi(end-1,2:end-1).*W(end,2:end-1).*h.*1) ...
            + sum(XSi(2:end-1,2)    .*U(2:end-1,1)  .*h.*1) ...
            - sum(XSi(2:end-1,end-1).*U(2:end-1,end).*h.*1);

dsumCFedt   = sum(CFe(2,2:end-1)    .*W(1,2:end-1)  .*h.*1) ...
            - sum(CFe(end-1,2:end-1).*W(end,2:end-1).*h.*1) ...
            + sum(CFe(2:end-1,2)    .*U(2:end-1,1)  .*h.*1) ...
            - sum(CFe(2:end-1,end-1).*U(2:end-1,end).*h.*1);

dsumCSidt   = sum(CSi(2,2:end-1)    .*W(1,2:end-1)  .*h.*1) ...
            - sum(CSi(end-1,2:end-1).*W(end,2:end-1).*h.*1) ...
            + sum(CSi(2:end-1,2)    .*U(2:end-1,1)  .*h.*1) ...
            - sum(CSi(2:end-1,end-1).*U(2:end-1,end).*h.*1);

if step>1
HST.dM(stp) = (alpha2*HST.dM(stp-1) + alpha3*HST.dM(max(1,stp-2)) + (beta1*dsumMdt + beta2*dsumMdto + beta3*dsumMdtoo)*dt)/alpha1; 
HST.dS(stp) = (alpha2*HST.dS(stp-1) + alpha3*HST.dS(max(1,stp-2)) + (beta1*dsumSdt + beta2*dsumSdto + beta3*dsumSdtoo)*dt)/alpha1; 
HST.dXFe(stp) = (alpha2*HST.dXFe(stp-1) + alpha3*HST.dXFe(max(1,stp-2)) + (beta1*dsumXFedt + beta2*dsumXFedto + beta3*dsumXFedtoo)*dt)/alpha1; 
HST.dXSi(stp) = (alpha2*HST.dXSi(stp-1) + alpha3*HST.dXSi(max(1,stp-2)) + (beta1*dsumXSidt + beta2*dsumXSidto + beta3*dsumXSidtoo)*dt)/alpha1; 
HST.dCFe(stp) = (alpha2*HST.dCFe(stp-1) + alpha3*HST.dCFe(max(1,stp-2)) + (beta1*dsumCFedt + beta2*dsumCFedto + beta3*dsumCFedtoo)*dt)/alpha1; 
HST.dCSi(stp) = (alpha2*HST.dCSi(stp-1) + alpha3*HST.dCSi(max(1,stp-2)) + (beta1*dsumCSidt + beta2*dsumCSidto + beta3*dsumCSidtoo)*dt)/alpha1; 
else
    HST.dM  (stp) = 0;
    HST.dS  (stp) = 0;
    HST.dXFe(stp) = 0;
    HST.dXSi(stp) = 0;
    HST.dCFe(stp) = 0;
    HST.dCSi(stp) = 0;
end

% Record conservation error (mass - mass change)/initial mass
HST.EM  (t)     = (HST.sumM  (t) - HST.dM  (t))./HST.sumM  (1) - 1;
HST.ES  (t)     = (HST.sumS  (t) - HST.dS  (t))./HST.sumS  (1) - 1;
HST.EXFe(t)     = (HST.sumXFe(t) - HST.dXFe(t))./HST.sumXFe(1) - 1;
HST.EXSi(t)     = (HST.sumXSi(t) - HST.dXSi(t))./HST.sumXSi(1) - 1;
HST.ECSi(t)     = (HST.sumCSi(t) - HST.dCSi(t))./HST.sumCSi(1) - 1;
HST.ECFe(t)     = (HST.sumCFe(t) - HST.dCFe(t))./HST.sumCFe(1) - 1;

if Nx <= 10 && Nz <= 10  % create 0D plots
HST.T(stp) = mean(mean(T));
HST.xFe(stp) = mean(mean(xFe));
HST.xSi(stp) = mean(mean(xSi));
HST.cFe(stp) = mean(mean(cFe));
HST.cSi(stp) = mean(mean(cSi));
HST.flFe(stp)= mean(mean(flFe));
HST.fsFe(stp)= mean(mean(fsFe));
HST.flSi(stp)= mean(mean(flSi));
HST.fsSi(stp)= mean(mean(fsSi));
HST.rho(stp)= mean(mean(rho));
HST.Eta(stp)= mean(mean(Eta));
end

