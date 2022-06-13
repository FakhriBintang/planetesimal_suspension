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

    % non-linear iteration loop
    while resnorm/resnorm0 >= NUM.reltol && resnorm >= NUM.abstol && iter <= NUM.maxit || iter < 2

        if NUM.step<=1; NUM.theta = 1; else; NUM.theta = 0.5; end

        % solve thermo-chemical equations
        solve_thermochem;
        % Temporary: Check latent and sensible heat proportions

%         if ~mod(NUM.step,RUN.nop)
%             masss = sum(sum(MAT.rho(2:end-1,2:end-1))) .*NUM.h.*NUM.h.*1;
%             Q_sen = PHY.Cp.*MAT.rho.*(SOL.T-To);
%             Q_lat = MAT.Ds.*(MAT.rho-rhoo);
%             Q = Q_sen + Q_lat;
%             TX = {'Interpreter','Latex'}; FS = {'FontSize',12};
%             TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',10};
%             UN = {'Units','Centimeters'};
% 
%             figure(30); clf
%             subplot(1,2,1)
%             plot(mean(Q_sen(2:end-1,2:end-1),2),NUM.zP(2:end-1),'LineWidth',2); axis ij tight; box on; hold on
%             plot(mean(Q_lat(2:end-1,2:end-1),2),NUM.zP(2:end-1),'LineWidth',2);
%             plot(mean(Q(2:end-1,2:end-1),2),NUM.zP(2:end-1),'LineWidth',2);
%             legend('$Q_{sens}$', '$Q_{lat}$', '$\bar{Q}$',TX{:},FS{:}); set(gca,TL{:},TS{:})
%             title('$Q_{sens}, Q_{lat}, \bar{Q}$',TX{:},FS{:}); set(gca,TL{:},TS{:});
%             ylabel('Depth [m]',TX{:},FS{:});
% 
%             subplot(1,2,2)
%             plot(mean(Q(2:end-1,2:end-1),2),NUM.zP(2:end-1),'LineWidth',2); axis ij tight; box on; hold on
%             plot(mean(SOL.H(2:end-1,2:end-1)-Ho(2:end-1,2:end-1),2),NUM.zP(2:end-1),'LineWidth',2);
%             legend('$Q_{sum}$', '$Q_{dH}$',TX{:},FS{:}); set(gca,TL{:},TS{:})
%             title('$Q_{sum}, Q_{dH}$',TX{:},FS{:}); set(gca,TL{:},TS{:});
%             ylabel('Depth [m]',TX{:},FS{:});
%         end

        if RUN.rad; radioactive_decay; end

        % update non-linear parameters and auxiliary variables
        up2date;


        if ~mod(iter,1) || iter ==0
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
%         figure(18)
%         subplot(1,3,1)
%         plot(mean(SOL.H(2:end-1,2:end-1),2),NUM.zP(2:end-1),'LineWidth',2); axis ij tight; box on;
%         title('$H [J]$',TX{:},FS{:}); set(gca,TL{:},TS{:});
%         ylabel('Depth [m]',TX{:},FS{:});
%         subplot(1,3,2)
%         plot(mean(MAT.Ds(2:end-1,2:end-1),2),NUM.zP(2:end-1),'LineWidth',2); axis ij tight; box on;
%         title('$Ds$',TX{:},FS{:}); set(gca,TL{:},TS{:});
%         ylabel('Depth [m]',TX{:},FS{:});
%         subplot(1,3,3)
%         plot(mean(SOL.H(2:end-1,2:end-1)-Ho(2:end-1,2:end-1),2),NUM.zP(2:end-1),'LineWidth',2); axis ij tight; box on;
%         title('$\Delta H$',TX{:},FS{:}); set(gca,TL{:},TS{:});
%         ylabel('Depth [m]',TX{:},FS{:});
%         %temp
        if RUN.rad
            figure(21); plot(HST.time/NUM.yr, HIST.Hr(2:end)); xlabel ('time'); ylabel('Hr (W/m^3)')
        end
    end

    % increment time
    NUM.step = NUM.step + 1;
    NUM.time = NUM.time + NUM.dt;

end