% planetesimal: main model routine

% print run header
fprintf(1,'\n\n************************************************************\n');
fprintf(1,    '*****  planetesimal  |  %s  |  %s  *****\n'         ,RUN.ID,datetime);
fprintf(1,    '************************************************************\n\n');


% %% initialise model run
% initialise;

S   = [SOL.W(:);SOL.U(:);SOL.P(:)];

while NUM.time <= NUM.tend && NUM.step <= NUM.maxstep
    % print time step header
    fprintf(1,'\n*****  step = %d;  dt = %1.4e;  time = %1.4e yr \n\n',NUM.step,NUM.dt/NUM.yr,NUM.time/NUM.yr);
    
    % store previous solutions
    MAT.Rhoo = MAT.Rhot;
    To                  = SOL.T;     % store previous temperature solution
    dTdto               = SOL.dTdt;  % store previous rate of change
    phio                = SOL.phi;
    
    % reset residuals and iteration count
    resnorm  = 1e3;
    resnorm0 = resnorm;
    iter     = 0;
    
        % non-linear iteration loop
    while resnorm/resnorm0 >= NUM.restol && resnorm >= NUM.abstol && iter <= NUM.maxit || iter <= 2
            
        % solve thermo-chemical equations
        solve_thermochem;
        
        % update non-linear parameters and auxiliary variables
        up2date;
        
        % solve fluid-mechanics equations
        solve_fluidmech;
        
        % update non-linear parameters and auxiliary variables
        up2date;
        
%         % report convergence
%         report;

        iter = iter+1;
    end
    
    if ~mod(NUM.step,round(2*RUN.nop/NUM.CFL))
        output;   
    end
    
    % increment time
    NUM.step = NUM.step + 1;
    NUM.time = NUM.time + NUM.dt;
    
end