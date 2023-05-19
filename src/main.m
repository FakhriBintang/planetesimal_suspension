% planetesimal: main model routine

% print run header
fprintf(1,'\n\n************************************************************\n');
fprintf(1,    '*****  planetesimal  |  %s  |  %s  *****\n'         ,RunID,datetime);
fprintf(1,    '************************************************************\n\n');


% %% initialise model run
initialise;

while time <= tend && step <= maxstep
    % print time step header
    fprintf(1,'\n*****  step = %d;  dt = %1.4e;  time = %1.4e yr;  %s\n\n',step,dt/yr,time/yr,dtlimit);

%     figure(100);clf

    if     strcmp(TINT,'bwei') || step==1 % first step / 1st-order backward-Euler implicit scheme
        alpha1 = 1; alpha2 = 1; alpha3 = 0;
        beta1  = 1; beta2  = 0; beta3  = 0;
    elseif strcmp(TINT,'bd3i') || step==2 % second step / 2nd-order 3-point backward-difference implicit scheme
        alpha1 = 3/2; alpha2 = 4/2; alpha3 = -1/2;
        beta1  = 1;   beta2  =  0;  beta3  = 0;
    elseif strcmp(TINT,'cnsi')            % other steps / 2nd-order Crank-Nicolson semi-implicit scheme
        alpha1 = 1;   alpha2 = 1;   alpha3 = 0;
        beta1  = 1/2; beta2  = 1/2; beta3  = 0;
    elseif strcmp(TINT,'bd3s')            % other steps / 2nd-order 3-point backward-difference semi-implicit scheme
        alpha1 = 3/2; alpha2 = 4/2; alpha3 = -1/2;
        beta1  = 3/4; beta2  = 2/4; beta3  = -1/4;
    end

    % store previous solution and auxiliary fields
    Soo     = So;       So    = S;
    XFeoo   = XFeo;     XFeo  = XFe;
    XSioo   = XSio;     XSio  = XSi;
    CSioo   = CSio;     CSio  = CSi;
    CFeoo   = CFeo;     CFeo  = CFe;
    FsFeoo   = FsFeo;   FsFeo  = FsFe;
    FsSioo   = FsSio;   FsSio  = FsSi;
    dSdtoo  = dSdto;    dSdto = dSdt;
    dXFedtoo= dXFedto;  dXFedto = dXFedt;
    dXSidtoo= dXSidto;  dXSidto = dXSidt;
    dCFedtoo= dCFedto;  dCFedto = dCFedt;
    dCSidtoo= dCSidto;  dCSidto = dCSidt;
    dFFedtoo= dFFedto;  dFFedto = dFFedt;
    dFSidtoo= dFSidto;  dFSidto = dFSidt;
    rhooo   = rhoo;     rhoo    = rho;
    Div_rhoVoo = Div_rhoVo; Div_rhoVo = Div_rhoV;
    Div_Vo  = Div_V;
    dto     = dt;
    % temp
    cFeo = cFe; cSio = cSi;

    % temp for radioactive decay
    NAlo  = NAl;
    dNdto = dNdt;

    % reset residuals and iteration count
    resnorm   = 1;
    resnorm0  = resnorm;
    iter      = 0;


    % non-linear iteration loop
    while resnorm/resnorm0 >= reltol && resnorm >= abstol && iter <= maxit

        % solve thermo-chemical equations
        solve_thermochem;

%         if radheat; radioactive_decay; else; Hr = Hr0; end

        % update non-linear parameters and auxiliary variables
        up2date;

        % solve fluid-mechanical equations
        solve_fluidmech;

        if ~bnchm
            resnorm = resnorm_TC + resnorm_VP;
            if iter == 0
                resnorm0 = resnorm+TINY;
            end
            fprintf(1,'  ---  it = %d;  abs res = %1.4e;  rel res = %1.4e  \n',iter,resnorm,resnorm/resnorm0)

%             figure(100); if iter==1; clf; else; hold on; end
%             plot(iter,log10(resnorm_TC),'b.',iter,log10(resnorm_VP),'r.',iter,log10(resnorm),'k.','MarkerSize',15,'LineWidth',1.5); box on; axis tight;
%             drawnow;
        end

        iter = iter+1;
    end

    % update mass and energy errors
    history;


    fprintf(1,'         min T   =  %4.1f;    mean T   = %4.1f;    max T   = %4.1f;   [degC]\n' ,min(T(:)),mean(T(:)-273.15),max(T(:)));
    fprintf(1,'         min x_{Fe}   =  %1.4f;    mean x_{Fe}   = %1.4f;    max x_{Fe}   = %1.4f;   [wt]\n'   ,min(xFe(:)  ),mean(xFe(:)  ),max(xFe(:)  ));
    fprintf(1,'         min c_{Fe}   =  %1.4f;    mean c_{Fe}   = %1.4f;    max c_{Fe}   = %1.4f;   [wt]\n'   ,min(cFe(:)  ),mean(cFe(:)  ),max(cFe(:)  ));
    fprintf(1,'         min c_{Si}   =  %1.4f;    mean c_{Si}   = %1.4f;    max c_{Si}   = %1.4f;   [wt]\n'   ,min(cSi(:)  ),mean(cSi(:)  ),max(cSi(:)  ));
    
    fprintf(1,'         min sumX=  %4.1f;    mean sumX= %4.1f;    max sumX= %4.1f;   [kg/m3]\n'  ,min(XFe(:)+XSi(:)),mean(XFe(:)+XSi(:))   ,max(XFe(:)+XSi(:)));
    fprintf(1,'         min rho =  %4.1f;    mean rho = %4.1f;    max rho = %4.1f;   [kg/m3]\n'  ,min(rho(:)),mean(rho(:))   ,max(rho(:)));
    if ~mod(step,nop) %round(2*nop/CFL))
        output;
        %        suppfigs; % supplementary figures for testing
%         if radheat
%             figure(21); plot(HST.time/yr, HIST.Hr(2:end)); xlabel ('time'); ylabel('Hr (W/m^3)')
%         end
    end

    % break script for NaN
    if isnan(resnorm)
        output
        print('Error, Breaking script')
        %         output
        return
    end
    % increment time
    dt;
    step = step + 1;
    time = time + dt;

end