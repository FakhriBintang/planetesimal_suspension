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
    Ho        = SOL.H;
    To        = SOL.T;
    fFelo     = CHM.fFel;
    fSilo     = CHM.fSil;
    cFeo      = CHM.cFe;
    cSio      = CHM.cSi;
    CSio      = CHM.CSi;
    CFeo      = CHM.CFe;
    XFeo      = CHM.XFe;
    xFeo      = CHM.xFe;
    xSio      = CHM.xSi;
    dHdto     = dHdt;
    dXdto     = dXdt;
    dCSidto   = dCSidt;
    dCFedto   = dCFedt;
    dfFedto   = dfFedt;
    dfSidto   = dfSidt;
    Pto       = SOL.Pt;

    % temp for radioactive decay
    NAlo = NAl;
    dNdto = dNdt;
    
    % reset residuals and iteration count
    resnorm   = 1e3;
    resnorm0  = resnorm;
    iter      = 0;
    tic
    % non-linear iteration loop
    while resnorm/resnorm0 >= NUM.reltol && resnorm >= NUM.abstol && iter <= NUM.maxit || iter < 2
        
        if NUM.step<=1; NUM.theta = 1; else; NUM.theta = 0.5; end
        
        % solve thermo-chemical equations
        solve_thermochem;
        if RUN.rad; radioactive_decay; end
     
        % update non-linear parameters and auxiliary variables
        up2date;
        

        if ~mod(iter,5) || iter ==0
        % solve fluid-mechanics equations
        solve_fluidmech;
        % update non-linear parameters and auxiliary variables
        up2date;
        end
        
        if ~RUN.bnchm
            resnorm = resnorm_TC + resnorm_VP;
            if iter == 0
                resnorm0 = resnorm+TINY;
            end
            fprintf(1,'  ---  it = %d;  abs res = %1.4e;  rel res = %1.4e  \n',iter,resnorm,resnorm/resnorm0)
            
           
            figure(100)
            plot(iter, log10(resnorm), 'k.','MarkerSize',20); axis xy tight; box on; hold on;
            drawnow;
            
        end
        
        iter = iter+1;
    end
toc
    % update mass and energy errors
    history;
    %temporary
    HIST.Hr = [HIST.Hr mean(MAT.Hr(:))];

    
    if ~mod(NUM.step,RUN.nop) %round(2*RUN.nop/NUM.CFL))
        output;   
            %temp
            if RUN.rad
    figure(21); plot(HST.time/NUM.yr, HIST.Hr(2:end)); xlabel ('time'); ylabel('Hr (W/m^3)')
            end
    end
    
    % increment time
    NUM.step = NUM.step + 1;
    NUM.time = NUM.time + NUM.dt;
    
end