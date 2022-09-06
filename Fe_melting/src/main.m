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
    So        = SOL.S;
    dSdto     = dSdt;
    dHdto     = dHdt;
    dXdto     = dXdt;
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

    % non-linear iteration loop
    while resnorm/resnorm0 >= NUM.reltol && resnorm >= NUM.abstol && iter <= NUM.maxit || iter < 2

        if NUM.step<=1; NUM.theta = 1; else; NUM.theta = 0.5; end

        % solve thermo-chemical equations
        if RK3
            solve_thermochem2;
        else
            solve_thermochem_entr;
%             solve_thermochem;
        end
        if RUN.rad; radioactive_decay; end

        % update non-linear parameters and auxiliary variables
        up2date;


        if ~mod(iter,1) || iter ==0
            % solve fluid-mechanics equations
%             solve_fluidmech_2;
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


            figure(100); if iter==1; clf; else; hold on; end
            plot(iter,log10(resnorm_TC),'b.',iter,log10(resnorm_VP),'r.',iter,log10(resnorm),'k.','MarkerSize',15,'LineWidth',1.5); box on; axis tight;
            drawnow;

        end

        iter = iter+1;
    end


    % calculate rayleigh number
    Ray = (NUM.D^3)*(max(MAT.rhoFes(:))-min(MAT.rho(:)))*0.1/mean(MAT.Eta(:))/(mean(mean(MAT.kT./MAT.rho./MAT.rho./PHY.Cp)));
    if Ray>1e7
        disp(['WARNING: Rayleigh Number too high, log10(Ra) = ', num2str(log10(Ray))])
    end

    % update mass and energy errors
    history;
    %temporary
    HIST.Hr = [HIST.Hr mean(MAT.Hr(:))];


    if ~mod(NUM.step,RUN.nop) %round(2*RUN.nop/NUM.CFL))
        output;
%         figure(32); clf
%         subplot(1,4,1)
%         plot(mean(SOL.S(2:end-1,2:end-1),2)  ,NUM.zP(2:end-1),'-k','LineWidth',2); axis ij tight; box on; hold on
%         plot(mean(sumS(2:end-1,2:end-1),2)   ,NUM.zP(2:end-1),'--k','LineWidth',2);
%         title('S')
%         subplot(1,4,2)
%         plot(mean(dSdt(2:end-1,2:end-1),2)  ,NUM.zP(2:end-1),'-k','LineWidth',2); axis ij tight; box on; hold on
%         plot(mean(advn_S(2:end-1,2:end-1),2),NUM.zP(2:end-1),'-r','LineWidth',2);
%         plot(mean(diff_T(2:end-1,2:end-1),2),NUM.zP(2:end-1),'-g','LineWidth',2);
%         plot(mean(diss_T(2:end-1,2:end-1),2),NUM.zP(2:end-1),'-b','LineWidth',2);
%         legend('dSdt', 'adv S', 'diff T', 'diss H')
%         subplot(1,4,3)
%         plot(mean(SOL.T(2:end-1,2:end-1),2)  ,NUM.zP(2:end-1),'-k','LineWidth',2); axis ij tight; box on; hold on
%         title('T')
%         subplot(1,4,4)
%         plot(mean(SOL.sFes(2:end-1,2:end-1),2),NUM.zP(2:end-1),'LineWidth',2); axis ij tight; box on; hold on
%         plot(mean(SOL.sFel(2:end-1,2:end-1),2),NUM.zP(2:end-1),'LineWidth',2);
%         plot(mean(SOL.sSis(2:end-1,2:end-1),2),NUM.zP(2:end-1),'LineWidth',2);
%         plot(mean(SOL.sSil(2:end-1,2:end-1),2),NUM.zP(2:end-1),'LineWidth',2);
        if RUN.rad
            figure(21); plot(HST.time/NUM.yr, HIST.Hr(2:end)); xlabel ('time'); ylabel('Hr (W/m^3)')
        end
    end

    % increment time
    NUM.dt
    NUM.step = NUM.step + 1;
    NUM.time = NUM.time + NUM.dt;

end