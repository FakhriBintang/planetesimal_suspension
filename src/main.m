% planetesimal: main model routine
while time <= tend && step <= maxstep
    % print time step header
    fprintf(1,'\n*****  step = %d;  dt = %1.4e;  time = %1.4e yr;  %s\n\n',step,dt/yr,time/yr,dtlimit);

%     figure(100);clf

    if     strcmp(TINT,'bwei') || step==1 % first step / 1st-order backward-Euler implicit scheme
        a1 = 1; a2 = 1; a3 = 0;
        b1  = 1; b2  = 0; b3  = 0;
    elseif strcmp(TINT,'bd3i') || step==2 % second step / 2nd-order 3-point backward-difference implicit scheme
        a1 = 3/2; a2 = 4/2; a3 = -1/2;
        b1  = 1;   b2  =  0;  b3  = 0;
    elseif strcmp(TINT,'cnsi')            % other steps / 2nd-order Crank-Nicolson semi-implicit scheme
        a1 = 1;   a2 = 1;   a3 = 0;
        b1  = 1/2; b2  = 1/2; b3  = 0;
    elseif strcmp(TINT,'bd3s')            % other steps / 2nd-order 3-point backward-difference semi-implicit scheme
        a1 = 3/2; a2 = 4/2; a3 = -1/2;
        b1  = 3/4; b2  = 2/4; b3  = -1/4;
    end

    % store previous solution and auxiliary fields
    if radheat
    n26Aloo     = n26Alo;   n26Alo      = n26Al; 
    dndtoo      = dndto;    dndto       = dndt; 
    end
    Soo         = So;       So          = S;
    XFeoo       = XFeo;     XFeo        = XFe;
    XSioo       = XSio;     XSio        = XSi;
    CSioo       = CSio;     CSio        = CSi;
    CFeoo       = CFeo;     CFeo        = CFe;
    FsFeoo      = FsFeo;    FsFeo       = FsFe;
    FsSioo      = FsSio;    FsSio       = FsSi;
    FlFeoo      = FlFeo;    FlFeo       = FlFe;
    FlSioo      = FlSio;    FlSio       = FlSi;
    dSdtoo      = dSdto;    dSdto       = dSdt;
    dXFedtoo    = dXFedto;  dXFedto     = dXFedt;
    dXSidtoo    = dXSidto;  dXSidto     = dXSidt;
    dCFedtoo    = dCFedto;  dCFedto     = dCFedt;
    dCSidtoo    = dCSidto;  dCSidto     = dCSidt;
    dFsFedtoo   = dFsFedto; dFsFedto    = dFsFedt;
    dFsSidtoo   = dFsSidto; dFsSidto    = dFsSidt;
    dFlFedtoo   = dFlFedto; dFlFedto    = dFlFedt;
    dFlSidtoo   = dFlSidto; dFlSidto    = dFlSidt;
    rhooo       = rhoo;     rhoo        = rho;
    Div_rhoVoo  = Div_rhoVo;Div_rhoVo   = Div_rhoV;
    Div_Vo      = Div_V;
    dto         = dt;
    % temp
    cFeo = cFe; cSio = cSi;

    % reset residuals and iteration count
    resnorm   = 1;
    resnorm0  = resnorm;
    iter      = 0;


    % non-linear iteration loop
    while resnorm/resnorm0 >= reltol && resnorm >= abstol && iter <= maxit
        % solve thermo-chemical equations
        solve_thermochem;

        % solve fluid-mechanical equations
        solve_fluidmech;

        % update non-linear parameters and auxiliary variables
        up2date;

        if ~bnchm
            resnorm = resnorm_TC + resnorm_VP;
            if iter == 0
                resnorm0 = resnorm+TINY;
            end
            fprintf(1,'  ---  it = %d;  abs res = %1.4e;  rel res = %1.4e  \n',iter,resnorm,resnorm/resnorm0)

            figure(100); if iter==1; clf; else; hold on; end
            % plot(iter,log10(resnorm_TC),'b.',iter,log10(resnorm_VP),'r.',iter,log10(resnorm),'k.','MarkerSize',15,'LineWidth',1.5); box on; axis tight;
            plot(iter,log10(resnorm_TC),'b.','MarkerSize',15,'LineWidth',1.5); box on; axis tight;
            drawnow;
        end
        iter = iter+1;
    end

    % update mass and energy errors
    history;
    Re     = D*rho(2:end-1,2:end-1).*Vel(2:end-1,2:end-1)./Eta(2:end-1,2:end-1)./10;


    fprintf(1,'         min T   =  %4.1f;    mean T   = %4.1f;    max T   = %4.1f;   [deg k]\n' ,min(T(:)),mean(T(:)),max(T(:)));
    fprintf(1,'         min x_{Fe}   =  %1.4f;    mean x_{Fe}   = %1.4f;    max x_{Fe}   = %1.4f;   [wt]\n'   ,min(xFe(:)  ),mean(xFe(:)  ),max(xFe(:)  ));
    fprintf(1,'         min c_{Fe}   =  %1.4f;    mean c_{Fe}   = %1.4f;    max c_{Fe}   = %1.4f;   [wt]\n'   ,min(cFe(:)  ),mean(cFe(:)  ),max(cFe(:)  ));
    fprintf(1,'         min c_{Si}   =  %1.4f;    mean c_{Si}   = %1.4f;    max c_{Si}   = %1.4f;   [wt]\n'   ,min(cSi(:)  ),mean(cSi(:)  ),max(cSi(:)  ));
    fprintf(1,'         min Re   =  %1.4f;    mean Re   = %1.4f;    max Re   = %1.4f;   [wt]\n'   ,min(Re(:)  ),mean(Re(:)  ),max(Re(:)  ));
    fprintf(1,'         min sumX=  %4.1f;    mean sumX= %4.1f;    max sumX= %4.1f;   [kg/m3]\n'  ,min(XFe(:)+XSi(:)),mean(XFe(:)+XSi(:))   ,max(XFe(:)+XSi(:)));
    fprintf(1,'         min rho =  %4.1f;    mean rho = %4.1f;    max rho = %4.1f;   [kg/m3]\n'  ,min(rho(:)),mean(rho(:))   ,max(rho(:)));
    if ~mod(step,nop) %round(2*nop/CFL))
        output
    end

    % break script for NaN
    if isnan(resnorm)
        output
        print('Error, Breaking script')
        return
    end
    % increment time
    dt;
    step = step + 1;
    time = time + dt;

end