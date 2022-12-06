% planetesimal: main model routine

% print run header
fprintf(1,'\n\n************************************************************\n');
fprintf(1,    '*****  planetesimal  |  %s  |  %s  *****\n'         ,RUN.ID,datetime);
fprintf(1,    '************************************************************\n\n');


% %% initialise model run
initialise;

% run manufactured solution benchmark on fluid mechanics solver if specified
if RUN.bnchm
    mms; return;
end

while NUM.time <= NUM.tend && NUM.step <= NUM.maxstep 
    % print time step header
    fprintf(1,'\n*****  step = %d;  dt = %1.4e;  time = %1.4e yr;  %s\n\n',NUM.step,NUM.dt/NUM.yr,NUM.time/NUM.yr,dtlimit);

    figure(100);clf

    % store previous solution and auxiliary fields
    rhoo      = MAT.rho;
    To        = SOL.T;
    fFelo     = CHM.fFel;
    fSilo     = CHM.fSil;
    cFeo      = CHM.cFe;
    cSio      = CHM.cSi;
    CSio      = CHM.CSi;
    CFeo      = CHM.CFe;
    XFeo      = CHM.XFe;
    XSio      = CHM.XSi;
    FFeo      = FFe;
    FSio      = FSi;
    xFeo      = CHM.xFe;
    xSio      = CHM.xSi;
    So        = SOL.S;
    dSdto     = dSdt;
    dXFedto   = dXFedt;
    dXSidto   = dXSidt;
    dCSidto   = dCSidt;
    dCFedto   = dCFedt;
    dFFedto   = dFFedt;
    dFSidto   = dFSidt;
    Pto       = SOL.Pt;

    % temp for radioactive decay
    NAlo = NAl;
    dNdto = dNdt;

    % reset residuals and iteration count
    resnorm   = 1e3;
    resnorm0  = resnorm;
    iter      = 0;

    if NUM.step<=1; NUM.theta = 1; else; NUM.theta = 0.5; end

    % non-linear iteration loop
    while resnorm/resnorm0 >= NUM.reltol && resnorm >= NUM.abstol && iter <= NUM.maxit || iter < 2

        % solve thermo-chemical equations
            solve_thermochem;

        if RUN.rad; radioactive_decay; end

        % update non-linear parameters and auxiliary variables
        up2date;


        if ~mod(iter,1) || iter ==0
            % solve fluid-mechanics equations
            solve_fluidmech;
%             % update non-linear parameters and auxiliary variables
%             up2date;
        end

        if ~RUN.bnchm
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

    if ~mod(NUM.step,RUN.nop) %round(2*RUN.nop/NUM.CFL))
        output;
%        suppfigs; % supplementary figures for testing
        if RUN.rad
            figure(21); plot(HST.time/NUM.yr, HIST.Hr(2:end)); xlabel ('time'); ylabel('Hr (W/m^3)')
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
    NUM.dt;
    NUM.step = NUM.step + 1;
    NUM.time = NUM.time + NUM.dt;

end