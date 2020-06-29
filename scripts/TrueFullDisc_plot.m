%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script to plot True Fully-Discrete Fourier Analysis data for DG schemes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

plt_dump_flag=0;

P=2;           % Polynomial order
Prk = 3;       % RK order
Beta = 1.00;      % upwind_parameter
CFL_max = 0.209;

ratio = 0.5;
CFL = ratio .* CFL_max;

Niter= [1,10,100];
Nt= length(Niter);
[Kt,S,G, Kp,Sp,Gp, ~, ~]=DG_TrueFullDiscFourierAnalys(P,Prk....
    ,Beta,CFL_max,ratio,Niter,0);

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
ccolor_map = [0,0,0; 0,0,1; 0,1,1; 1,0,1; 0,1,0; 1,0,0;  0.5 1.0 0.5]; 
line_style={':','-.','-','--','-.','-','--',':'};
marker_type= {'o','>','*','s','<','+','v'};
%line_style={'--','-','-.',':','-.','-','--',':'}; % for beta=0, DGp2
% marker_type= {'o','>','s','*','<','+','v'};

if(P==5)
    marker_step=150;
    marker_indices = [1,2,5:marker_step:length(Kp)-marker_step,length(Kp)];
else
   marker_step=8;
   marker_indices = [5,marker_step:marker_step:length(Kp)-marker_step,length(Kp)-marker_step*0.5,length(Kp)]; 
end

kk = pi/4;
index(1) = find(Kp./(P+1) <= kk , 1, 'last' );
index(2) = find(Kt./(P+1) <= kk , 1, 'last' );
index

legend_fontsize=21;
axes_labels_fontsize=24;
axes_ticks_fontsize=22;
%=========================================================================%
%                  Dissipationn compare physical and true
%=========================================================================%
legend_temp = cell(1,2*Nt);
h=figure;
ax = gca;
ax.TickLabelInterpreter = 'latex';

for i=1:Nt
    pl(i)=plot(Kt(2:end)./N,G(i,2:end)....
        ,'color',ccolor_map(i,:),'linewidth',2....
        ,'LineStyle',char(line_style(i))); 
    hold on
    legend_temp(1,i) = { strcat('n= ',num2str(Niter(i)) )};
end

for i=1:Nt
    pl(i+Nt)=plot(Kp./N,Gp(i,:)...
        ,'color',ccolor_map(i,:),'linewidth',2....
        ,'LineStyle',char(line_style(i))....
        ,'Marker',char(marker_type(i))....
        ,'MarkerIndices',marker_indices,'MarkerSize',8); 
    hold on
    legend_temp(1,i+Nt) = { strcat('n= ',num2str(Niter(i)) )};
end
% pl(Nt+1) = plot(Kp./N, ones(1,length(Kp)),':k','linewidth',1.2);
legend_temp(Nt+1) = {'exact'};

ax.LineWidth = 1.0;
ax.GridLineStyle = '-';
% ax.GridColor = 'k';
ax.GridAlpha = 0.2;

xx=xlabel('$$K$$','Interpreter','latex','FontSize',axes_labels_fontsize);
yy=ylabel('$$G(K,\tau_{p,n})$$','Interpreter','latex','FontSize',axes_labels_fontsize);
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
set(gca,'Color',[0.9 0.9 0.9]);
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 7 5];

outerpos = ax.OuterPosition;
ti = ax.TightInset; 
ti(4)=0.03;
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
    err = abs(1-G(i,:));
    pl(i)=semilogy(Kt(2:end)./N,err(2:end)....
        ,'color',ccolor_map(i,:),'linewidth',2....
        ,'LineStyle',char(line_style(i))); 
    hold on
    legend_temp(1,i) = { strcat('n= ',num2str(Niter(i)) )};
end

for i=1:Nt
    err_p = abs(1-Gp(i,:));
    pl(i+Nt)=semilogy(Kp(3:end)./N,err_p(3:end)...
        ,'color',ccolor_map(i,:),'linewidth',2....
        ,'LineStyle',char(line_style(i))....
        ,'Marker',char(marker_type(i))....
        ,'MarkerIndices',marker_indices,'MarkerSize',8); 
    hold on
    legend_temp(1,i+Nt) = { strcat('\phi_{mode1}, n= ',num2str(Niter(i)) )};
end

ax.LineWidth = 1.0;
ax.GridLineStyle = '-';
% ax.GridColor = 'k';
ax.GridAlpha = 0.2;

xx=xlabel('$$K$$','Interpreter','latex','FontSize',axes_labels_fontsize);
yy=ylabel('$$|1-G(K,\tau_{p,n})|$$'....
    ,'Interpreter','latex','FontSize',axes_labels_fontsize);
legend(pl(1:Nt),legend_temp(1:Nt),'FontSize',legend_fontsize....
    ,'Location','southeast','Interpreter','latex')
grid on
ax.YMinorGrid = 'off';
xlim([0,pi/2])
ylim([1e-8,1])
xticks([0,pi/4,pi/2,3*pi/4,pi])
xticklabels({'\fontsize{21} 0','\fontsize{30}\pi\fontsize{21}/4'....
,'\fontsize{30}\pi\fontsize{21}/2'...
,'\fontsize{21}3\fontsize{30}\pi\fontsize{21}/4'...
,'\fontsize{30}\pi'})
yticks([1e-10,1e-8,1e-6,1e-4,1e-2,1,1e2,1e4])
ax.FontSize=axes_ticks_fontsize;
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
    err = abs(S(i,:))./(P+1);
    pl(i)=semilogy(Kt(2:end)./N,err(2:end)'....
        ,'color',ccolor_map(i,:),'linewidth',2....
        ,'LineStyle',char(line_style(i))); 
    hold on
    legend_temp(1,i) = { strcat('n= ',num2str(Niter(i)) )};
end

for i=1:Nt
    err1 = abs(Kp.*Niter(i).*CFL-Sp(i,:)')./(P+1);
    %err1 = abs(Kp-(Sp(i,:)'./CFL))./(P+1);
    pl(i+Nt)=semilogy(Kp./N,err1...
        ,'color',ccolor_map(i,:),'linewidth',2....
        ,'LineStyle',char(line_style(i))....
        ,'Marker',char(marker_type(i))....
        ,'MarkerIndices',marker_indices,'MarkerSize',8);
    hold on
    legend_temp(1,i+Nt) = { strcat('\phi_{mode1}, n= ',num2str(Niter(i)) )};
end

ax.LineWidth = 1.0;
ax.GridLineStyle = '-';
% ax.GridColor = 'k';
ax.GridAlpha = 0.2;

xx=xlabel('$$K$$','Interpreter','latex','FontSize',axes_labels_fontsize);
yy=ylabel('$$\Delta \psi(K,\tau_{p,n})$$'....
    ,'Interpreter','latex','FontSize',axes_labels_fontsize);
legend(pl(1:Nt),legend_temp(1:Nt),'FontSize',legend_fontsize....
    ,'Location','southeast','Interpreter','latex')
grid on
ax.YMinorGrid = 'off';
xlim([0,pi])
xticks([0,pi/4,pi/2,3*pi/4,pi])
xticklabels({'\fontsize{21} 0','\fontsize{30}\pi\fontsize{21}/4'....
,'\fontsize{30}\pi\fontsize{21}/2'...
,'\fontsize{21}3\fontsize{30}\pi\fontsize{21}/4'...
,'\fontsize{30}\pi'})
ylim([1e-10,100])
yticks([1e-10,1e-8,1e-6,1e-4,1e-2,1,1e2,1e4])
ax.FontSize=axes_ticks_fontsize;
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
