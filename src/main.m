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

    figure(100);clf

    % store previous solution and auxiliary fields
    So        = S;
    XFeo      = XFe;
    XSio      = XSi;
    CSio      = CSi;
    CFeo      = CFe;
    FFeo      = FsFe;
    FSio      = FsSi;
    dSdto     = dSdt;
    dXFedto   = dXFedt;
    dXSidto   = dXSidt;
    dCSidto   = dCSidt;
    dCFedto   = dCFedt;
    dFFedto   = dFFedt;
    dFSidto   = dFSidt;
    rhoo      = rho;
    Div_rhoVo = Div_rhoV;
    % temp
    cFeo = cFe; cSio = cSi;

    % temp for radioactive decay
    NAlo  = NAl;
    dNdto = dNdt;

    % reset residuals and iteration count
    resnorm   = 1e3;
    resnorm0  = resnorm;
    iter      = 0;


    % non-linear iteration loop
    while resnorm/resnorm0 >= reltol && resnorm >= abstol && iter <= maxit || iter < 2

        % solve thermo-chemical equations
        solve_thermochem;

        if radheat; radioactive_decay; else; Hr = Hr0; end

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

            figure(100); if iter==1; clf; else; hold on; end
            plot(iter,log10(resnorm_TC),'b.',iter,log10(resnorm_VP),'r.',iter,log10(resnorm),'k.','MarkerSize',15,'LineWidth',1.5); box on; axis tight;
            drawnow;
        end

        iter = iter+1;
    end

    % update mass and energy errors
    history;

    if ~mod(step,nop) %round(2*nop/CFL))
        output;
%        suppfigs; % supplementary figures for testing
        if radheat
            figure(21); plot(HST.time/yr, HIST.Hr(2:end)); xlabel ('time'); ylabel('Hr (W/m^3)')
        end
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