% t = step;

dsumMdtoo  = dsumMdto; dsumMdto = dsumMdt;
dsumSdtoo = dsumSdto; dsumSdto   = dsumSdt;
dsumXFedtoo = dsumXFedto; dsumXFedto = dsumXFedt;
dsumXSidtoo = dsumXSidto; dsumXSidto = dsumXSidt;
dsumCFedtoo = dsumCFedto; dsumCFedto = dsumCFedt;
dsumCSidtoo = dsumCSidto; dsumCSidto = dsumCSidt;

stp = max(1,step);

HST.time(stp) = time;

% record conserved masses at current timestep
if radheat
HST.n26Al (stp)      = n26Al; end

% m0      = rho0*4*pi()/3*rw(end)^3; % mass of bottom boundary
HST.sumM  (stp)      = sum(RHO(2:end-1,2)*4*pi()/3.*(rw(1:end-1).^3-rw(2:end).^3)); % mass of sphere above boundary
HST.sumS  (stp)      = sum((  S(2:end-1,2)+S0(2:end-1,2))...
                                         *4*pi()/3.*(rw(1:end-1).^3-rw(2:end).^3));
HST.sumXFe(stp)      = sum(XFe(2:end-1,2)*4*pi()/3.*(rw(1:end-1).^3-rw(2:end).^3));
HST.sumXSi(stp)      = sum(XSi(2:end-1,2)*4*pi()/3.*(rw(1:end-1).^3-rw(2:end).^3));
HST.sumCFe(stp)      = sum(XFe(2:end-1,2)*4*pi()/3.*(rw(1:end-1).^3-rw(2:end).^3));
HST.sumCSi(stp)      = sum(XSi(2:end-1,2)*4*pi()/3.*(rw(1:end-1).^3-rw(2:end).^3));


% record expected rates of volume change imposed by variable boundaries
% calculated relative to surface area of the outer boundary sphere
SArea       = 4*pi()*rw(1)^2;
dsumMdt     = rho(2,2)*WBG(1).*SArea;
dsumSdt     = (S(2,2) + Hr(2,2) + diss_T(2:2)).*WBG(1).*SArea;
dsumXFedt   = XFe(2,2)*WBG(1).*SArea;
dsumXSidt   = XSi(2,2)*WBG(1).*SArea;
dsumCFedt   = CFe(2,2)*WBG(1).*SArea;
dsumCSidt   = CSi(2,2)*WBG(1).*SArea;

% if stp>=2; HST.DM(stp  ) = HST.DM(max(1,stp-1)  ) + dsumMdt*dt; else; HST.DM(stp  ) = 0; end  % [kg]
% if stp>=2; HST.DS(stp  ) = HST.DS(max(1,stp-1)  ) + dsumSdt*dt; else; HST.DS(stp  ) = 0; end  % [kg]
if step>=2
HST.dM  (stp) = (a2*HST.dM  (stp-1) + a3*HST.dM  (max(1,stp-2)) + (b1*dsumMdt   + b2*dsumMdto   + b3*dsumMdtoo  )*dt)/a1; 
HST.dS  (stp) = (a2*HST.dS  (stp-1) + a3*HST.dS  (max(1,stp-2)) + (b1*dsumSdt   + b2*dsumSdto   + b3*dsumSdtoo  )*dt)/a1; 
HST.dXFe(stp) = (a2*HST.dXFe(stp-1) + a3*HST.dXFe(max(1,stp-2)) + (b1*dsumXFedt + b2*dsumXFedto + b3*dsumXFedtoo)*dt)/a1; 
HST.dXSi(stp) = (a2*HST.dXSi(stp-1) + a3*HST.dXSi(max(1,stp-2)) + (b1*dsumXSidt + b2*dsumXSidto + b3*dsumXSidtoo)*dt)/a1; 
HST.dCFe(stp) = (a2*HST.dCFe(stp-1) + a3*HST.dCFe(max(1,stp-2)) + (b1*dsumCFedt + b2*dsumCFedto + b3*dsumCFedtoo)*dt)/a1; 
HST.dCSi(stp) = (a2*HST.dCSi(stp-1) + a3*HST.dCSi(max(1,stp-2)) + (b1*dsumCSidt + b2*dsumCSidto + b3*dsumCSidtoo)*dt)/a1; 
else
    HST.dM  (stp) = 0;
    HST.dS  (stp) = 0;
    HST.dXFe(stp) = 0;
    HST.dXSi(stp) = 0;
    HST.dCFe(stp) = 0;
    HST.dCSi(stp) = 0;
end

% Record conservation error (mass - mass change)/initial mass
HST.EM  (stp)     = (HST.sumM  (stp) - HST.dM  (stp))./HST.sumM  (1) - 1;
HST.ES  (stp)     = (HST.sumS  (stp) - HST.dS  (stp))./HST.sumS  (1) - 1;
HST.EXFe(stp)     = (HST.sumXFe(stp) - HST.dXFe(stp))./HST.sumXFe(1) - 1;
HST.EXSi(stp)     = (HST.sumXSi(stp) - HST.dXSi(stp))./HST.sumXSi(1) - 1;
HST.ECSi(stp)     = (HST.sumCSi(stp) - HST.dCSi(stp))./HST.sumCSi(1) - 1;
HST.ECFe(stp)     = (HST.sumCFe(stp) - HST.dCFe(stp))./HST.sumCFe(1) - 1;

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

