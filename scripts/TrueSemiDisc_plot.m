%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script to plot True Semi-Discrete Fourier Analysis data for DG schemes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

plt_dump_flag=1;

P=5;           % Polynomial order
Beta = 1.00;      % upwind_parameter

% t= [0.1,1,10];
tau_p=[0.1,0.5,6];
taup_name={'0.1','0.5','6.00'};
t=tau_p./(P+1);
Nt= length(t);
parfor ii=1:length(t)
    [Kt(ii,:),S(ii,:),G(ii,:), Kp(ii,:),Sp(ii,:),Gp(ii,:), ~,~]=....
        DG_TrueSemiDiscFourierAnalys(P,Beta,t(ii),0);
end

N=P+1;

%Preparing directories:
%------------------------
%Preparing directories:
%------------------------
resdir='../results/';
outdir=strcat(resdir,'DGp',num2str(P),'/');
datadir=strcat(outdir,'data/');
figdir=strcat(outdir,'figures/');
png_fig = strcat(figdir,'png/');
eps_fig = strcat(figdir,'eps/');
if(~isdir(resdir))
    command=strcat('mkdir',{' '},resdir);
    system(command{1});
end
if(~isdir(outdir))
    command=strcat('mkdir',{' '},outdir);
    system(command{1});
end
if(~isdir(datadir))
    command=strcat('mkdir',{' '},datadir);
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

ccolor_map = hsv(length(Nt));
ccolor_map = [0,0,0; 0,0,1; 0,1,1; 0,1,0; 1,0,0; 1,0,1; 0.5 1.0 0.5 ; 1,0.7,1; 0,0.2,0.2; 0,0.5,0;]; 
line_style={'-.','-','--',':','-.','-','--',':','-','-.','--'};
marker_type= {'o','>','s','*','<','+','v','d','x','^'};
marker_step=10;
marker_indices = [5:marker_step:length(Kp)];
if(P==5)
    marker_step=20;
    marker_indices = [1,2,5:marker_step:length(Kp)];
end
legend_fontsize=21;
axes_labels_fontsize=24;
axes_ticks_fontsize=22;

%=========================================================================%
%                  Dissipationn compare physical and true
%=========================================================================%
legend_temp = cell(1,Nt);
h=figure;
ax = gca;
ax.TickLabelInterpreter = 'latex';

% pl(1) = plot(Kp(1,:)./N, ones(1,length(Kp(1,:))),':k','linewidth',1.0); hold on

for i=1:Nt
    pl(i)=plot(Kt(i,2:end)./N,G(i,2:end)....
        ,'color',ccolor_map(i,:),'linewidth',2....
        ,'LineStyle',char(line_style(i))); hold on
    legend_temp(1,i) = { strcat('$\tau_{p}= ',taup_name{i},'$')};
end

for i=1:Nt
    pl(i+Nt)=plot(Kp(i,:)./N,Gp(i,:)...
        ,'color',ccolor_map(i,:),'linewidth',2....
        ,'LineStyle',char(line_style(i))....
        ,'Marker',char(marker_type(i))....
        ,'MarkerIndices',marker_indices,'MarkerSize',8);
    hold on
%     legend_temp(1,i+Nt) = { strcat('$\tau_{p}= ',num2str(tau_p(i)),'$')};
end
% pl(Nt+1) = plot(Kp(i,:)./N, ones(1,length(Kp)),':k','linewidth',0.5);
%legend_temp(Nt+1) = {'Exact, G=1.0'};

ax.LineWidth = 1.0;
ax.GridLineStyle = '-';
% ax.GridColor = 'k';
ax.GridAlpha = 0.2;

xx=xlabel('$$K$$','Interpreter','latex','FontSize',axes_labels_fontsize);
yy=ylabel('$$G(K,\tau_{p})$$','Interpreter','latex','FontSize',axes_labels_fontsize);
legend(pl(1:Nt),legend_temp(1:Nt),'FontSize',legend_fontsize....
    ,'Location','southwest','Interpreter','latex')
% legend([pl(1),pl(3)],{'True, combined-mode','Physical-mode'}....
%     ,'FontSize',legend_fontsize,'Location','southwest'....
%     ,'Interpreter','latex')
grid on
xlim([0,pi])
ylim([-0.05,1.1])
xticks([0,pi/4,pi/2,3*pi/4,pi])
xticklabels({'\fontsize{21} 0','\fontsize{30}\pi\fontsize{21}/4'....
,'\fontsize{30}\pi\fontsize{21}/2'...
,'\fontsize{21}3\fontsize{30}\pi\fontsize{21}/4'...
,'\fontsize{30}\pi'})
ax.FontSize=axes_ticks_fontsize;
% ax.YAxis.FontSize=14;
set(gca,'Color',[0.9 0.9 0.9]);
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 7 5];

outerpos = ax.OuterPosition;
ti = ax.TightInset;
ti(4)=0.01;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

if(plt_dump_flag==1)
    ffname = strcat('p',num2str(P),'_beta',num2str(Beta),'_sdisc_G');
    fname_png = strcat(outdir,png_fig,ffname);
    fname_eps = strcat(outdir,eps_fig,ffname);
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
    err = abs(1-G(i,:));
    pl(i)=semilogy(Kt(i,1:end)./N,err(1:end),'color',ccolor_map(i,:)...
        ,'linewidth',2,'LineStyle',char(line_style(i))); hold on
%     legend_temp(1,i) = { strcat('$\tau_{p}= ',num2str(tau_p(i)),'$')};
end

for i=1:Nt
    err_p = abs(1-Gp(i,:));
    pl(i+Nt)=semilogy(Kp(i,1:end)./N,err_p(1:end)...
        ,'color',ccolor_map(i,:),'linewidth',2....
        ,'LineStyle',char(line_style(i))....
        ,'Marker',char(marker_type(i))....
        ,'MarkerIndices',marker_indices,'MarkerSize',8); hold on
%     legend_temp(1,i+Nt) = { strcat('$\phi_{mode1}, \tau_{p}= ',num2str(tau_p(i)),'$')};
end

ax.LineWidth = 1.0;
ax.GridLineStyle = '-';
% ax.GridColor = 'k';
ax.GridAlpha = 0.2;

xx=xlabel('$$K$$','Interpreter','latex','FontSize',axes_labels_fontsize);
yy=ylabel('$$|1-G(K,\tau_{p})|$$'....
    ,'Interpreter','latex','FontSize',axes_labels_fontsize);
legend(pl(1:Nt),legend_temp(1:Nt),'FontSize',legend_fontsize....
    ,'Location','southeast','Interpreter','latex')
grid on
xlim([0,pi/2])
ylim([1e-8,1])
xticks([0,pi/4,pi/2,3*pi/4,pi])
xticklabels({'\fontsize{21} 0','\fontsize{30}\pi\fontsize{21}/4'....
,'\fontsize{30}\pi\fontsize{21}/2'...
,'\fontsize{21}3\fontsize{30}\pi\fontsize{21}/4'...
,'\fontsize{30}\pi'})
yticks([1e-10,1e-8,1e-6,1e-4,1e-2,1,1e2,1e4])
ax.FontSize=axes_ticks_fontsize;
% ax.YAxis.FontSize=14;
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
    ffname = strcat('p',num2str(P),'_beta',num2str(Beta),'_sdisc_G_err');
    fname_png = strcat(outdir,png_fig,ffname);
    fname_eps = strcat(outdir,eps_fig,ffname);
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
    err = abs(S(i,:));
    pl(i)=semilogy(Kt(i,2:end)./N,err(2:end)','color',ccolor_map(i,:)...
        ,'linewidth',1.8,'LineStyle',char(line_style(i))); hold on
%     legend_temp(1,i) = { strcat('$\tau_{p}= ',num2str(tau_p(i)),'$')};
end

for i=1:Nt
    err1 = abs(Kp(i,:).*t(i)-Sp(i,:));
    pl(i+Nt)=semilogy(Kp(i,:)./N,err1...
        ,'color',ccolor_map(i,:),'linewidth',1.8....
        ,'LineStyle',char(line_style(i))....
        ,'Marker',char(marker_type(i))....
        ,'MarkerIndices',marker_indices,'MarkerSize',8); hold on
%     legend_temp(1,i+Nt) = { strcat('$\phi_{mode1}, \tau_{p}= ',num2str(tau_p(i)),'$')};
end

ax.LineWidth = 1.0;
ax.GridLineStyle = '-';
% ax.GridColor = 'k';
ax.GridAlpha = 0.2;

xx=xlabel('$$K$$','Interpreter','latex','FontSize',axes_labels_fontsize);
yy=ylabel('$$\Delta \psi(K,\tau_{p})$$'....
    ,'Interpreter','latex','FontSize',axes_labels_fontsize);
legend(pl(1:Nt),legend_temp(1:Nt),'FontSize',legend_fontsize....
    ,'Location','southeast','Interpreter','latex')
grid on
xlim([0,pi])
xticks([0,pi/4,pi/2,3*pi/4,pi])
xticklabels({'\fontsize{21} 0','\fontsize{30}\pi\fontsize{21}/4'....
,'\fontsize{30}\pi\fontsize{21}/2'...
,'\fontsize{21}3\fontsize{30}\pi\fontsize{21}/4'...
,'\fontsize{30}\pi'})
ylim([1e-10,100])
yticks([1e-10,1e-8,1e-6,1e-4,1e-2,1,1e2,1e4])
ax.FontSize=axes_ticks_fontsize;
% ax.YAxis.FontSize=14;
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
    ffname = strcat('p',num2str(P),'_beta',num2str(Beta),'_sdisc_S');
    fname_png = strcat(outdir,png_fig,ffname);
    fname_eps = strcat(outdir,eps_fig,ffname);
    print(fname_png,'-dpng','-r0')
    print(fname_eps,'-depsc','-r0')
end
