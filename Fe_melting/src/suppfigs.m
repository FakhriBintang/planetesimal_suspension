% Supplementary figures for error testing

fh32 = figure(32); clf
        subplot(1,4,1)
        plot(mean(SOL.S(2:end-1,2:end-1),2)  ,NUM.zP(2:end-1),'-k','LineWidth',2); axis ij tight; box on; hold on
        plot(mean(sumS(2:end-1,2:end-1),2)   ,NUM.zP(2:end-1),'--k','LineWidth',2);
        title('S')
        subplot(1,4,2)
        plot(mean(dSdt,2)  ,NUM.zP(2:end-1),'-k','LineWidth',2); axis ij tight; box on; hold on
        plot(mean(advn_S,2),NUM.zP(2:end-1),'-r','LineWidth',2);
        plot(mean(diff_S,2),NUM.zP(2:end-1),'-g','LineWidth',2);
        plot(mean(diss_T,2),NUM.zP(2:end-1),'-b','LineWidth',2);
        legend('dSdt', 'adv S', 'diff T', 'diss H')
        subplot(1,4,3)
        plot(mean(SOL.T(2:end-1,2:end-1),2)  ,NUM.zP(2:end-1),'-k','LineWidth',2); axis ij tight; box on; hold on
        title('T')
        subplot(1,4,4)
        plot(mean(SOL.sFes(2:end-1,2:end-1),2),NUM.zP(2:end-1),'LineWidth',2); axis ij tight; box on; hold on
        plot(mean(SOL.sFel(2:end-1,2:end-1),2),NUM.zP(2:end-1),'LineWidth',2);
        plot(mean(SOL.sSis(2:end-1,2:end-1),2),NUM.zP(2:end-1),'LineWidth',2);
        plot(mean(SOL.sSil(2:end-1,2:end-1),2),NUM.zP(2:end-1),'LineWidth',2);

        fh33 = figure(33); clf
        subplot(1,7,1)
        plot(mean(CHM.cFe(2:end-1,2:end-1),2) - mean(cFeo(2:end-1,2:end-1),2)  ,NUM.zP(2:end-1),'-r','LineWidth',2); axis ij tight; box on; hold on
        plot(mean(CHM.cSi(2:end-1,2:end-1),2) - mean(cSio(2:end-1,2:end-1),2)  ,NUM.zP(2:end-1),'-b','LineWidth',2);
        legend('\Delta cFe', '\Delta cSi')
        title('dc')
        subplot(1,7,2)
        plot(mean(CHM.CFe(2:end-1,2:end-1),2) - mean(CFeo(2:end-1,2:end-1),2)  ,NUM.zP(2:end-1),'-r','LineWidth',2); axis ij tight; box on; hold on
        plot(mean(CHM.CSi(2:end-1,2:end-1),2) - mean(CSio(2:end-1,2:end-1),2)  ,NUM.zP(2:end-1),'-b','LineWidth',2);
        legend('\Delta CFe', '\Delta CSi')
        title('dC')
        subplot(1,7,3)
        plot(mean(dCFedt,2) ,NUM.zP(2:end-1),'-r','LineWidth',2); axis ij tight; box on; hold on
        plot(mean(dCSidt,2) ,NUM.zP(2:end-1),'-b','LineWidth',2);
        legend('dCFedt', 'dCSidt')
        title('dCdt')

        subplot(1,7,4)
        plot(mean(CFeo(2:end-1,2:end-1),2)  ,NUM.zP(2:end-1),'--r','LineWidth',1); axis ij tight; box on; hold on
        plot(mean(CHM.CFe(2:end-1,2:end-1),2)  ,NUM.zP(2:end-1),'-r','LineWidth',2);
        legend('CFeo', 'CFe')
        subplot(1,7,5)
        plot(mean(cFeo(2:end-1,2:end-1),2)  ,NUM.zP(2:end-1),'--r','LineWidth',1); axis ij tight; box on; hold on
        plot(mean(CHM.cFe(2:end-1,2:end-1),2)  ,NUM.zP(2:end-1),'-r','LineWidth',2);
        legend('cFeo', 'cFe')

        subplot(1,7,6)
        plot(mean(CSio(2:end-1,2:end-1),2)  ,NUM.zP(2:end-1),'--b','LineWidth',1); axis ij tight; box on; hold on
        plot(mean(CHM.CSi(2:end-1,2:end-1),2)  ,NUM.zP(2:end-1),'-b','LineWidth',2);
        legend('CSio', 'CSi')
        subplot(1,7,7)
        plot(mean(cSio(2:end-1,2:end-1),2)  ,NUM.zP(2:end-1),'--b','LineWidth',1); axis ij tight; box on; hold on
        plot(mean(CHM.cSi(2:end-1,2:end-1),2)  ,NUM.zP(2:end-1),'-b','LineWidth',2);
        legend('cSio', 'cSi')

fh34 = figure(34); clf
%         subplot(1,4,1)
%         plot(mean(dFFedt(2:end-1,2:end-1),2)  ,NUM.zP(2:end-1),'-k','LineWidth',2); axis ij tight; box on; hold on
%         plot(mean(CHM.GFe(2:end-1,2:end-1),2)   ,NUM.zP(2:end-1),'-b','LineWidth',1.5);
% %         plot(mean(advn_FFe(2:end-1,2:end-1),2)   ,NUM.zP(2:end-1),'-r','LineWidth',1.5);
%         title('dFFedt')
%         legend('dFFedt', '$\Gamma_{Fe}$','adv ')
%         subplot(1,4,2)
%         plot(mean(dFSidt(2:end-1,2:end-1),2)  ,NUM.zP(2:end-1),'-k','LineWidth',2); axis ij tight; box on; hold on
%         plot(mean(CHM.GSi(2:end-1,2:end-1),2)   ,NUM.zP(2:end-1),'-b','LineWidth',1.5);
% %         plot(mean(advn_FSi(2:end-1,2:end-1),2)   ,NUM.zP(2:end-1),'-r','LineWidth',1.5);
%         title('dFSidt')
%         legend('dFSidt', '$\Gamma_{Si}$','adv ')
%         if ~XSolve
        subplot(1,5,1)
        plot(mean(CHM.XSi(2:end-1,2:end-1),2)  ,NUM.zP(2:end-1),'-k','LineWidth',2); axis ij tight; box on; hold on
        title('XSi')
        subplot(1,5,2)
        plot(mean(CHM.XFe(2:end-1,2:end-1),2)  ,NUM.zP(2:end-1),'-k','LineWidth',2); axis ij tight; box on; hold on
        title('XFe')
        subplot(1,5,3)
        plot(mean(MAT.rho(2:end-1,2:end-1),2) - mean(CHM.XSi(2:end-1,2:end-1),2) - mean(CHM.XFe(2:end-1,2:end-1),2)  ,NUM.zP(2:end-1),'-k','LineWidth',2); axis ij tight; box on; hold on
        title('Xerror')
        subplot(1,5,4)
        plot(1 - mean(CHM.xSi(2:end-1,2:end-1),2) - mean(CHM.xFe(2:end-1,2:end-1),2)  ,NUM.zP(2:end-1),'-k','LineWidth',2); axis ij tight; box on; hold on;
        title('xerror')

        subplot(1,5,5)
        plot(mean(dXFedt,2)  ,NUM.zP(2:end-1),'-r','LineWidth',2); axis ij tight; box on; hold on;
%         if XSolve
%            plot(mean(dXSidt(2:end-1,2:end-1),2)  ,NUM.zP(2:end-1),'-b','LineWidth',2);
%         else
%             % reconstruct Chemical density rate of change
%             dXSidt = (CHM.XSi - XSio)./NUM.dt;
%             plot(mean(dXSidt(2:end-1,2:end-1),2)  ,NUM.zP(2:end-1),'-b','LineWidth',2);

%         end
%         legend('dXFedt', 'dXSidt')
%         title('dXdt')

fh39 = figure(39); clf
        subplot(1,4,1)
        plot(mean(dSdt,2)  ,NUM.zP(2:end-1),'-k','LineWidth',2); axis ij tight; box on; hold on
        title('dSdt')
        subplot(1,4,2)
        plot(mean(dXFedt,2)  ,NUM.zP(2:end-1),'-k','LineWidth',2); axis ij tight; box on; hold on
        title('dFedt')
        subplot(1,4,3)
        plot(mean(dCFedt,2)  ,NUM.zP(2:end-1),'-k','LineWidth',2); axis ij tight; box on; hold on
        plot(mean(dCSidt,2)  ,NUM.zP(2:end-1),'-k','LineWidth',2);
        title('dCFedt','dCSidt')
        subplot(1,4,4)
        plot(mean(dFFedt,2)  ,NUM.zP(2:end-1),'-k','LineWidth',2); axis ij tight; box on; hold on
        title('dFdt')


        if RUN.save

        name = [outpath '/',RUN.ID,'_dSdt_',num2str(floor(NUM.step/RUN.nop))]; % figure 4
        print(fh32,name,'-dpng','-r300','-opengl');

        name = [outpath '/',RUN.ID,'_dCdt_',num2str(floor(NUM.step/RUN.nop))]; % figure 4
        print(fh33,name,'-dpng','-r300','-opengl');

%         name = [outpath '/',RUN.ID,'_dXdt_',num2str(floor(NUM.step/RUN.nop))]; % figure 4
%         print(fh34,name,'-dpng','-r300','-opengl');
        end


%         subplot(1,5,5)
%         plot(1-mean(CHM.xSi(2:end-1,2:end-1),2) - mean(CHM.xFe(2:end-1,2:end-1),2)  ,NUM.zP(2:end-1),'-k','LineWidth',2); axis ij tight; box on; hold on


%         subplot(1,4,3)
%         plot(mean(SOL.T(2:end-1,2:end-1),2)  ,NUM.zP(2:end-1),'-k','LineWidth',2); axis ij tight; box on; hold on
%         title('T')
%         subplot(1,4,4)
%         plot(mean(SOL.sFes(2:end-1,2:end-1),2),NUM.zP(2:end-1),'LineWidth',2); axis ij tight; box on; hold on
%         plot(mean(SOL.sFel(2:end-1,2:end-1),2),NUM.zP(2:end-1),'LineWidth',2);
%         plot(mean(SOL.sSis(2:end-1,2:end-1),2),NUM.zP(2:end-1),'LineWidth',2);
%         plot(mean(SOL.sSil(2:end-1,2:end-1),2),NUM.zP(2:end-1),'LineWidth',2);