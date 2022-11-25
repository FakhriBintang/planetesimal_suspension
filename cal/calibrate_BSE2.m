clear all; close all
addpath('../src'); close all;

% calibration run options
runID     = 'Phsdg3';     % run ID for output files; [system name_wt.% SiO2_wt.% H2O]
holdfig   = 0;                           % set to 1 to hold figures, to 0 for new figures
linestyle = '-';                         % set line style for plots
save_plot = 0;                           % turn on (1) to save output file in /out directory

%% load data
% Load MELTS tables in csv format
Liquidcompupper = load('./Javoy2010/Liquid_comp_upper.csv');   % output file of MELTS for upper phase loop
Solidcompupper  = load('./Javoy2010/Solid_comp_upper.csv');   % output file of MELTS for upper phase loop
Bulkcompupper   = load('./Javoy2010/Bulk_comp_upper.csv');   % output file of MELTS for upper phase loop

Liquidcomplower = load('./McDonough_1995_1/Liquid_comp_lower.csv');   % output file of MELTS for upper phase loop
Solidcomplower  = load('./McDonough_1995_1/Solid_comp_lower.csv');   % output file of MELTS for upper phase loop
Bulkcomplower   = load('./McDonough_1995_1/Bulk_comp_lower.csv');   % output file of MELTS for upper phase loop
%% extract parameters
% upper phase loop
P_upper = Bulkcompupper(:,1).'.*1e5;
T_upper = Bulkcompupper(:,2).';
c_upper = Bulkcompupper(:,4).'./100;
% read in phase fractions and compositions from MELTS tables
m_upper  = Liquidcompupper(:,3).'./Bulkcompupper(:,3).';
x_upper  = Solidcompupper(:,3).'./Bulkcompupper(:,3).';
cm_upper = Liquidcompupper(:,4).'./100;
cx_upper = Solidcompupper(:,4).'./100;

% lower phase loop
P_lower = Bulkcomplower(:,1).'.*1e5;
T_lower = Bulkcomplower(:,2).';
c_lower = Bulkcomplower(:,4).'./100;
% read in phase fractions and compositions from MELTS tables
m_lower  = Liquidcomplower(:,3).'./Bulkcomplower(:,3).';
x_lower  = Solidcomplower(:,3).'./Bulkcomplower(:,3).';
cm_lower = Liquidcomplower(:,4).'./100;
cx_lower = Solidcomplower(:,4).'./100;

%% initial phase diagram parameters
cphs0_bst  =  0.4128;               % phase diagram lower bound composition [wt SiO2]
cphs1_bst  =  0.784;               % phase diagram upper bound composition [wt SiO2]
Tphs0_bst  =  980;               % phase diagram lower bound temperature [degC]
Tphs1_bst  =  1750;                % phase diagram upper bound temperature [degC]
PhDg_bst   =  [20,10,1.0,0.93];  % Phase diagram curvature factor (> 1)
perCm_bst  =  0.557;               % peritectic liquidus composition [wt SiO2]
perCx_bst  =  0.493;               % peritectic solidus  composition [wt SiO2]
perT_bst   =  1210;                % peritectic temperature [degC]
clap       =  0;                % Clapeyron slope for P-dependence of melting T [degC/Pa]
beta       =  0.9;                 % iterative lag parameter phase diagram [1]

%% run calibration loop
misfit  = 1e3; tol = 0.4;
bestfit = misfit;

time0 = tic;
timelimit = 60*60*0.25; % Limit 15min %%change back to 5min after completion of misfit

while misfit > tol
    % set phase diagram parameters
    cphs0    =  cphs0_bst * (1 + randn(1)*1e-3);          % phase diagram lower bound composition [wt SiO2]
    cphs1    =  cphs1_bst * (1 + randn(1)*1e-3);          % phase diagram upper bound composition [wt SiO2]
    Tphs0    =  Tphs0_bst * (1 + randn(1)*1e-3);          % phase diagram lower bound temperature [degC]
    Tphs1    =  Tphs1_bst * (1 + randn(1)*1e-3);          % phase diagram upper bound temperature [degC]
    PhDg     =  PhDg_bst .* (1 + randn(1,4)*1e-3);        % Phase diagram curvature factor (> 1)
    perCm    =  perCm_bst * (1 + randn(1)*3e-3);          % peritectic liquidus composition [wt SiO2]
    perCx    =  perCx_bst * (1 + randn(1)*3e-3);          % peritectic solidus  composition [wt SiO2]
    perT     =  perT_bst  * (1 + randn(1)*3e-3);          % peritectic temperature [degC]

    % equilibrium phase fractions and compositions - upperphsloop
    [xq_upper,cxq_upper,cmq_upper]  =  equilibrium(T_upper,c_upper,zeros(size(T_upper)),Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,PhDg,1e-16);
    mq_upper = 1-xq_upper;

    % equilibrium phase fractions and compositions - lowerphsloop
    [xq_lower,cxq_lower,cmq_lower]  =  equilibrium(T_lower,c_lower,zeros(size(T_lower)),Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,PhDg,1e-16);
    mq_lower = 1-xq_lower;


    misfit =          norm( m_upper.*linspace(1.15,1,length(m_upper))- mq_upper,2)./norm( mq_upper,2) ...
        + norm( x_upper- xq_upper,2)./norm( xq_upper,2) ...
        + norm(cm_upper-cmq_upper,2)./norm(cmq_upper,2) ...
        + norm(cx_upper-cxq_upper,2)./norm(cxq_upper,2)*0.8;
    misfit = misfit + norm( m_lower- mq_lower,2)./norm( mq_lower,2) ...
        + norm( x_lower- xq_lower,2)./norm( xq_lower,2) ...
        + norm(cm_lower-cmq_lower,2)./norm(cmq_lower,2) ...
        + norm(cx_lower-cxq_lower,2)./norm(cxq_lower,2)*0.8;


    if misfit < 1.0005*bestfit         % if misfit < bestfit set bestfit values as misfit
        fprintf(1,'   misfit = %1.4e \n',misfit);
        bestfit = misfit;       %set bestfit values of parameter to model parameters
        cphs0_bst = cphs0;
        cphs1_bst = cphs1;
        Tphs0_bst = Tphs0;
        Tphs1_bst = Tphs1;
        PhDg_bst  = PhDg;
        perCm_bst = perCm;
        perCx_bst = perCx;
        perT_bst  = perT;

        %        plot phase fractions
        figure(1); if ~holdfig; clf; end
        plot(T_upper,xq_upper.*100,'k-' ,T_upper,mq_upper.*100,'r-' , 'LineWidth',2); hold on; box on; axis tight;
        plot(T_lower,xq_lower.*100,'k-' ,T_lower,mq_lower.*100,'r-' , 'LineWidth',2);
        plot(T_upper,x_upper*100,'kd',T_upper,m_upper*100,'rd','LineWidth',2,'MarkerSize',5);
        plot(T_lower,x_lower*100,'kv',T_lower,m_lower*100,'rv','LineWidth',2,'MarkerSize',5);
        set(gca,'TickLabelInterpreter','latex','FontSize',13)
        title('Melting model','Interpreter','latex','FontSize',18)
        xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
        ylabel('Phase fractions [wt\%]','Interpreter','latex','FontSize',15)

        % plot major phase compositions
        figure(2); if ~holdfig; clf; end
        plot(T_upper,cxq_upper.*100,'b-' ,T_upper,cmq_upper.*100,'r-' ,'LineWidth',2); hold on; box on; axis tight;
        plot(T_lower,cxq_lower.*100,'b-' ,T_lower,cmq_lower.*100,'r-' ,'LineWidth',2);
        plot(T_upper,cx_upper.*100,'bd' ,T_upper,cm_upper.*100,'rd' ,'LineWidth',2);
        plot(T_lower,cx_lower.*100,'bv' ,T_lower,cm_lower.*100,'rv' ,'LineWidth',2);
        set(gca,'TickLabelInterpreter','latex','FontSize',13)
        title('Phase compositions','Interpreter','latex','FontSize',18)
        xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
        ylabel('Major component [wt\% SiO$_2$]','Interpreter','latex','FontSize',15)

        drawnow
    end

    if toc(time0)>timelimit  % if loop takes too long, break after time limit
        break
    end
end

cphs0 = cphs0_bst
cphs1 = cphs1_bst
Tphs0 = Tphs0_bst
Tphs1 = Tphs1_bst
PhDg  = PhDg_bst
perCm = perCm_bst
perCx = perCx_bst
perT  = perT_bst

% set ranges for control variables T, c, v, P - BSE
T = linspace(min(T_upper)-100,max(T_upper)+100,1e3);  % temperature range [degC]
c = 0.514 * ones(size(T));                       % major component range [wt SiO2]
v = 0.000 * ones(size(T));                       % volatile component range [wt H2O]
P = 0 * ones(size(T));                       % pressure range [Pa]

% equilibrium phase fractions and compositions - BSE
[xq,cxq,cmq]  =  equilibrium(T,c,zeros(size(T)),Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,PhDg,1e-16);


mq = 1-xq;

% plot phase diagram
figure(3); if ~holdfig; clf; end
vv = (4.8e-5.*mean(P(:)).^0.6 + 1e-9.*mean(P(:)))./100;
TT = linspace(Tphs0+mean(P(:))*clap,Tphs1+mean(P(:))*clap,1e3);
cc = [linspace(cphs1,(perCx+perCm)/2,ceil((perT-Tphs0)./(Tphs1-Tphs0)*1e3)),linspace((perCx+perCm)/2,cphs0,floor((perT-Tphs1)./(Tphs0-Tphs1)*1e3))];
[~,CCx,CCm]   = equilibrium(TT,cc,0.*TT,Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,PhDg,1e-16);
plot(CCx.*100,TT,'k-','LineWidth',2); axis tight; hold on; box on;
plot(CCm.*100,TT,'k-','LineWidth',2);
Tphs0s = Tphs0;
Tphs1s = Tphs1;
perTs  = perT;
TT = linspace(Tphs0s+mean(P(:))*clap,Tphs1s+mean(P(:))*clap,1e3);
cc = [linspace(cphs1,(perCx+perCm)/2,round((perTs-Tphs0s)./(Tphs1s-Tphs0s)*1e3)),linspace((perCx+perCm)/2,cphs0,round((perTs-Tphs1s)./(Tphs0s-Tphs1s)*1e3))];
[~,CCx,CCm]   = equilibrium(TT,cc,0.*TT,Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,PhDg,1e-16);
plot(CCx.*100,TT,'k-','LineWidth',2); axis tight; hold on; box on;
plot(CCm.*100,TT,'k-','LineWidth',2);
plot(cxq.*100,T-(P-mean(P(:))).*clap,'b','LineStyle',linestyle,'LineWidth',2);
plot(cmq.*100,T-(P-mean(P(:))).*clap,'r','LineStyle',linestyle,'LineWidth',2);
% plot(c./(1-fq+1e-16).*100,T-(P-mean(P(:))).*clap,'Color',[0.5 0.5 0.5],'LineStyle',linestyle,'LineWidth',2);
set(gca,'TickLabelInterpreter','latex','FontSize',13)
title('Phase Diagram','Interpreter','latex','FontSize',18)
xlabel('Major component [wt\% SiO$_2$]','Interpreter','latex','FontSize',15)
ylabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)




% create output directory
if ~isfolder(['../out/',runID])
    mkdir(['../out/',runID]);
end

% save output to file
if save_plot
    name = ['../out/',runID,'/',runID,'_phase_dgrm'];
    print(figure(1),name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_melt_model'];
    print(figure(2),name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_maj_compnt'];
    print(figure(3),name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_vol_compnt'];


    % save ('name','cphs0','cphs1','Tphs0','Tphs1','PhDg','perCm','perCx','perT','clap','dTH2O' 'beta', 'bestfit', '-ascii');
end

% for i = bestfit         %not yet done to generate an output file
%     X = []
%    fid = sprintf('output_matlab_%d.txt', i);
%
%     fprintf(fid,'%4.4f\n',[cphs0 cphs1 Tphs0 Tphs1 PhDg perCm perCx perT]);
%     fclose(fid, 'all');
% end