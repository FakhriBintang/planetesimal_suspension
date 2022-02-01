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
    
    figure(100);clf
    
    % store previous solutions
    MAT.rhoo = MAT.rhot;
    Kvbo = Kvb;
    Ho                  = SOL.H;     % store previous temperature solution
    To                  = SOL.T;
    cFeo                = CHM.cFe;
    cSio                = CHM.cSi;
    CSio                = CHM.CSi;
    CFeo                = CHM.CFe;
    dHdto               = dHdt;  % store previous rate of change
    XFeo                = CHM.XFe;
    dXdto               = dXdt;
    dCSidto             = dCSidt;
    dCFedto             = dCFedt;
    Pto                 = SOL.Pt;
    
%     phio                = SOL.phi;
    Div_rhovo           = Div_rhov;
    % reset residuals and iteration count
    resnorm  = 1e3;
    resnorm0 = resnorm;
    iter     = 0;
    
        % non-linear iteration loop
    while resnorm/resnorm0 >= NUM.restol && resnorm >= NUM.abstol && iter <= NUM.maxit || iter <= 2
%             for i = 1:10
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
% output
        if iter==0; resnorm0 = resnorm; end
        iter = iter+1;
%             end
    end

    % update mass error
    MassErr = MassErr - dVoldt.*NUM.dt;

    
    if ~mod(NUM.step,RUN.nop) %round(2*RUN.nop/NUM.CFL))
        output;   
    end
    
    % increment time
    NUM.step = NUM.step + 1;
    NUM.time = NUM.time + NUM.dt;
    
end