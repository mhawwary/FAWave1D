%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script to plot True Fourier Analysis data for multible schemes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc

plt_dump_flag=0;

P=5;           % Polynomial order
Prk = 4;       % RK order
Beta = [0.00,0.20,1.00];      % upwind_parameter
Nb = length(Beta);
CFL_max = [0.103, 0.108, 0.073];

ratio = 0.9;
CFL = ratio .* CFL_max;
Niter= [10];
Nt= length(Niter);
[Kt,S1,G1, Kp1,Sp1,Gp1, E1, Eex1]=DG_TrueFullDiscFourierAnalys(P,Prk....
    ,Beta(1),CFL_max(1),ratio,Niter,0);
[~,S2,G2, Kp2,Sp2,Gp2, E2, Eex2]=DG_TrueFullDiscFourierAnalys(P,Prk....
    ,Beta(2),CFL_max(2),ratio,Niter,0);
[~,S3,G3, Kp3,Sp3,Gp3, E3, Eex3]=DG_TrueFullDiscFourierAnalys(P,Prk....
    ,Beta(3),CFL_max(3),ratio,Niter,0);

N=P+1;

%Preparing directories:
%------------------------
outdir='../results/';
figdir=strcat(outdir,'figures/');
png_fig = strcat(figdir,'png/');
eps_fig = strcat(figdir,'eps/');
if(~isdir(outdir))
    command=strcat('mkdir',{' '},outdir);
    system(command{1});
end
if(~isdir(figdir))
    command=strcat('mkdir',{' '},figdir);
    system(command{1});
end
if(~isdir(png_fig))
    command=strcat('mkdir',{' '},png_fig);
    system(command{1});
end
if(~isdir(eps_fig))
    command=strcat('mkdir',{' '},eps_fig);
    system(command{1});
end
addpath(genpath('../src/'))

%=========================================================================%
%               Exact and True Energy as a function of time and K
%=========================================================================%
ccolor_map = hsv(length(Beta));
ccolor_map = [0,0,0; 0,0,1; 0,1,1; 0,1,0; 1,0,0; 1,0,1; 0.5 1.0 0.5]; 
% legend_temp = cell(1,2*Nt+1);
% h=figure;
% ax = gca;
% ax.TickLabelInterpreter = 'latex';
% 
% plot(Kt./N,Eex1(1,:),'color',ccolor_map(1,:)...
%         ,'linewidth',1.5,'LineStyle','--'), hold on
% legend_temp(1,1) = {'n=0'};
%     
% for i=1:Nt
%     plot(Kt./N,E1(i,:),'color',ccolor_map(i+1,:)...
%         ,'linewidth',1.5,'LineStyle','-'), hold on
%     legend_temp(1,i+1) = { strcat('n= ',num2str(Niter(i)) )};
% end
% 
% xx=xlabel('K/P+1','Interpreter','latex','FontSize',15);
% yy=ylabel('E'....
%     ,'Interpreter','latex','FontSize',15);
% legend(legend_temp(1:Nt+1),'FontSize',15....
%     ,'Location','southwest','Interpreter','latex')
% % title(strcat('FullyDiscrete, DGp2, "---" E$_{true}$, "- -" E$_{ex}$')....
% %     ,'FontSize',14,'FontWeight','Normal','Interpreter','latex')
% grid on
% xlim([0,pi])
% % ylim([0.4,1.1])
% xticks([0,pi/4,pi/2,3*pi/4,pi])
% xticklabels({'\fontsize{14} 0','\fontsize{20}\pi\fontsize{14}/4'....
% ,'\fontsize{20}\pi\fontsize{14}/2'...
% ,'\fontsize{14}3\fontsize{20}\pi\fontsize{14}/4'...
% ,'\fontsize{20}\pi'})
% ax.YAxis.FontSize=14;
% set(gca,'Color',[0.9 0.9 0.9]);
% h.PaperUnits = 'inches';
% h.PaperPosition = [0 0 7 5];
% 
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% 
% if(plt_dump_flag==1)
%     ffname = strcat('p',num2str(P),'_RK',num2str(Prk).....
%         ,'_beta',num2str(Beta),'_0',num2str(ratio*10),'CFLmax_E');
%     fname_png = strcat(png_fig,ffname);
%     fname_eps = strcat(eps_fig,ffname);
%     print(fname_png,'-dpng','-r0')
%     print(fname_eps,'-depsc','-r0')
% end
%=========================================================================%
%                  Dissipationn compare physical and true
%=========================================================================%
legend_temp = cell(1,2*Nb);
h=figure;
ax = gca;
ax.TickLabelInterpreter = 'latex';

% for i=1:Nt
%     pl(i)=plot(Kt./N,G1(i,:),'color',ccolor_map(i,:)...
%         ,'linewidth',1.5,'LineStyle','-'); hold on
%     legend_temp(1,i) = { strcat('n= ',num2str(Niter(i)) )};
% end
% 
% for i=1:Nt
%     pl(i+Nt)=plot(Kp1./N,Gp1(i,:)...
%     ,'color',ccolor_map(i,:),'linewidth',1.5,'LineStyle','--'); hold on
%     legend_temp(1,i+Nt) = { strcat('n= ',num2str(Niter(i)) )};
% end

i=1;
pl(i)=plot(Kt./N,G1,'color',ccolor_map(i,:)...
        ,'linewidth',1.5,'LineStyle','-'); hold on
legend_temp(1,i) = { strcat('$\beta= ',num2str(Beta(i)),'$')};
i=2;
pl(i)=plot(Kt./N,G2,'color',ccolor_map(i,:)...
        ,'linewidth',1.5,'LineStyle','-'); hold on
legend_temp(1,i) = { strcat('$\beta= ',num2str(Beta(i)),'$')};
i=3;
pl(i)=plot(Kt./N,G3,'color',ccolor_map(i,:)...
        ,'linewidth',1.5,'LineStyle','-'); hold on
legend_temp(1,i) = { strcat('$\beta= ',num2str(Beta(i)),'$')};
i=4;
pl(i)=plot(Kp1./N,Gp1...
    ,'color',ccolor_map(i-3,:),'linewidth',1.5,'LineStyle','--'); hold on
i=5;
pl(i)=plot(Kp2./N,Gp2...
    ,'color',ccolor_map(i-3,:),'linewidth',1.5,'LineStyle','--'); hold on
i=6;
pl(i)=plot(Kp3./N,Gp3...
    ,'color',ccolor_map(i-3,:),'linewidth',1.5,'LineStyle','--'); hold on


xx=xlabel('K/P+1','Interpreter','latex','FontSize',15);
yy=ylabel('$G^{n}$'....
    ,'Interpreter','latex','FontSize',15);
legend(pl(1:Nb),legend_temp(1:Nb),'FontSize',15....
    ,'Location','southwest','Interpreter','latex')
grid on
xlim([0,pi])
ylim([-0.05,1.1])
xticks([0,pi/4,pi/2,3*pi/4,pi])
xticklabels({'\fontsize{14} 0','\fontsize{20}\pi\fontsize{14}/4'....
,'\fontsize{20}\pi\fontsize{14}/2'...
,'\fontsize{14}3\fontsize{20}\pi\fontsize{14}/4'...
,'\fontsize{20}\pi'})
ax.YAxis.FontSize=14;
set(gca,'Color',[0.9 0.9 0.9]);
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 7 5];

outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

if(plt_dump_flag==1)
    ffname = strcat('p',num2str(P),'_RK',num2str(Prk).....
        ,'_beta',num2str(Beta),'_0',num2str(ratio*10),'CFLmax_G');
    fname_png = strcat(png_fig,ffname);
    fname_eps = strcat(eps_fig,ffname);
    print(fname_png,'-dpng','-r0')
    print(fname_eps,'-depsc','-r0')
end

%=========================================================================%
%                  Compare Dissipation error for physical and true 
%=========================================================================%
h=figure;
ax = gca;
ax.TickLabelInterpreter = 'latex';

for i=1:Nt
    err = abs(1-G1(i,:));
    pl(i)=semilogy(Kt(1:end)./N,err(1:end),'color',ccolor_map(i,:)...
        ,'linewidth',1.5,'LineStyle','-'); hold on
    legend_temp(1,i) = { strcat('n= ',num2str(Niter(i)) )};
end

for i=1:Nt
    err_p = abs(1-Gp1(i,:));
    pl(i+Nt)=semilogy(Kp1(3:end)./N,err_p(3:end)...
        ,'color',ccolor_map(i,:),'linewidth',1.5,'LineStyle','--'); hold on
    legend_temp(1,i+Nt) = { strcat('\phi_{mode1}, n= ',num2str(Niter(i)) )};
end

xx=xlabel('K/P+1','Interpreter','latex','FontSize',15);
yy=ylabel('$|1-G^{n}|$'....
    ,'Interpreter','latex','FontSize',15);
legend(pl(1:Nt),legend_temp(1:Nt),'FontSize',15....
    ,'Location','southeast','Interpreter','latex')
grid on
xlim([0,pi])
ylim([1e-10,100])
xticks([0,pi/4,pi/2,3*pi/4,pi])
xticklabels({'\fontsize{14} 0','\fontsize{20}\pi\fontsize{14}/4'....
,'\fontsize{20}\pi\fontsize{14}/2'...
,'\fontsize{14}3\fontsize{20}\pi\fontsize{14}/4'...
,'\fontsize{20}\pi'})
yticks([1e-10,1e-8,1e-6,1e-4,1e-2,1,1e2,1e4])
ax.YAxis.FontSize=14;
set(gca,'Color',[0.9 0.9 0.9]);
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 7 5];

outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

if(plt_dump_flag==1)
    ffname = strcat('p',num2str(P),'_RK',num2str(Prk).....
        ,'_beta',num2str(Beta),'_0',num2str(ratio*10),'CFLmax_G_err');
    fname_png = strcat(png_fig,ffname);
    fname_eps = strcat(eps_fig,ffname);
    print(fname_png,'-dpng','-r0')
    print(fname_eps,'-depsc','-r0')
end


%=========================================================================%
%                  Compare Dispersion error for physical and true 
%=========================================================================%
h=figure;
ax = gca;
ax.TickLabelInterpreter = 'latex';

for i=1:Nt
    err = abs(S1(i,:))./(P+1);
    pl(i)=semilogy(Kt./N,err,'color',ccolor_map(i,:)...
        ,'linewidth',1.5,'LineStyle','-'); hold on
    legend_temp(1,i) = { strcat('n= ',num2str(Niter(i)) )};
end

for i=1:Nt
    err1 = abs(Kp1.*Niter(i)*CFL-Sp1)./(P+1);
    pl(i+Nt)=semilogy(Kp1./N,err1...
        ,'color',ccolor_map(i,:),'linewidth',1.5,'LineStyle','--'); hold on
    legend_temp(1,i+Nt) = { strcat('\phi_{mode1}, n= ',num2str(Niter(i)) )};
end

xx=xlabel('K/P+1','Interpreter','latex','FontSize',15);
yy=ylabel('$\Delta \psi(t_{n})$'....
    ,'Interpreter','latex','FontSize',15);
legend(pl(1:Nt),legend_temp(1:Nt),'FontSize',15....
    ,'Location','southeast','Interpreter','latex')
grid on
xlim([0,pi])
xticks([0,pi/4,pi/2,3*pi/4,pi])
xticklabels({'\fontsize{14} 0','\fontsize{20}\pi\fontsize{14}/4'....
,'\fontsize{20}\pi\fontsize{14}/2'...
,'\fontsize{14}3\fontsize{20}\pi\fontsize{14}/4'...
,'\fontsize{20}\pi'})
ylim([1e-10,100])
yticks([1e-10,1e-8,1e-6,1e-4,1e-2,1,1e2,1e4])
ax.YAxis.FontSize=14;
set(gca,'Color',[0.9 0.9 0.9]);
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 7 5];

outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

if(plt_dump_flag==1)
    ffname = strcat('p',num2str(P),'_RK',num2str(Prk).....
        ,'_beta',num2str(Beta),'_0',num2str(ratio*10),'CFLmax_S');
    fname_png = strcat(png_fig,ffname);
    fname_eps = strcat(eps_fig,ffname);
    print(fname_png,'-dpng','-r0')
    print(fname_eps,'-depsc','-r0')
end
