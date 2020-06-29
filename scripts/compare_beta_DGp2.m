% Comparison of different DGp2 Beta Schemes:
close all
clear
clc

P=2;           % Polynomial order
Prk = 0;       % RK order
plt_dump_flag=0;

ccolor_map = [0,0,0; 0,0,1; 0,1,1; 1,0,1; 1,0,0; 0,1,0; 0.5 1.0 0.5]; 
%ccolor_map = [0,0,0; 0,0,0; 0,0,0; 0,0,0];  % black&white
line_style={'-','--',':','-.','-','-','--',':'};
marker_type= {'o','^','s','.','<','+','v'};
legend_fontsize=17;
axes_labels_fontsize=16;
axes_ticks_fontsize=16;


%Preparing directories:
%------------------------
resdir='../results/';
outdir = strcat(resdir,'DGp',num2str(P),'/');
datadir = strcat(outdir,'data/');
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

legend_wp={'$\beta=1.0$','$\beta=0.55$','$\beta=0.2$','$\beta=0.0$'};
legend_wd={'$\beta=1.0$','$\beta=0.55$','$\beta=0.2$'};

% Loading data:
Beta = 1.00;      % fully upwind flux
fname_wd = strcat(datadir,'DGp',num2str(P),'_Beta'...
    ,num2str(Beta),'_sdisc_wd.dat');
fname_wp = strcat(datadir,'DGp',num2str(P),'_Beta'...
    ,num2str(Beta),'_sdisc_wp.dat');


[K,wd,wp]=load_dominant_mode(P,0,Beta,0);
[K(:,2),wd(:,2),wp(:,2)]=load_dominant_mode(P,0,0.55,0);
[K(:,3),wd(:,3),wp(:,3)]=load_dominant_mode(P,0,0.2,0);
[K(:,4),wd(:,4),wp(:,4)]=load_dominant_mode(P,0,0.0,0);

marker_step = [11,8,6,5];
%marker_indices = 1:marker_step:length(K(:,1));

%------------------------
% Dissipation plot:
%------------------------
h=figure;
ax = gca;
ax.TickLabelInterpreter = 'latex';
for i=1:3
    marker_indices=[120:marker_step(i):length(K(:,i))-10];
    plot(K(:,i)./(P+1),wd(:,i)./(P+1),'color',ccolor_map(i,:),'linewidth',0.8....
                ,'LineStyle',char(line_style(i))....
                ,'Marker',char(marker_type(i))....
                ,'MarkerIndices',marker_indices,'MarkerSize',8....
                ),hold on
end
plot(K(:,1)./(P+1),zeros(1,length(K(:,1))),':k')

xlabel('$K$','Interpreter','latex','FontSize',14);
ylabel('$\mathcal{I}$m$(K_{m})$'...
    ,'Interpreter','latex','FontSize',14)
ylim([-4,1])
xlim([0,pi])
xticks([0,pi/4,pi/2,3*pi/4,pi])
xticklabels({'\fontsize{14} 0','\fontsize{20}\pi\fontsize{14}/4'....
,'\fontsize{20}\pi\fontsize{14}/2'...
,'\fontsize{14}3\fontsize{20}\pi\fontsize{14}/4'...
,'\fontsize{20}\pi'})
ax.YAxis.FontSize=14;

legend({'$\beta=1.00$','$\beta=0.55$','$\beta=0.20$','exact'}...
    ,'FontSize',15,'Location','southwest','Interpreter','latex')

h.PaperUnits = 'inches';
h.PaperPosition = [0 0 6 4];
set(gca,'Color',[0.9 0.9 0.9]);

outerpos = ax.OuterPosition;
ti = ax.TightInset; 
%ti(4) = 0.01;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

if(plt_dump_flag==1)
    ffname = strcat('p',num2str(P),'_compare_beta_sdisc_wd');
    fname_png = strcat(png_fig,ffname,'.png');
    fname_eps = strcat(eps_fig,ffname,'.eps');
    print(fname_png,'-dpng','-r0')
    print(fname_eps,'-depsc','-r0')
end


% Dispersion plot:
h=figure;
ax = gca;
ax.TickLabelInterpreter = 'latex';
for i=1:4
    marker_indices=[120:marker_step(i):length(K(:,i))-10];
    plot(K(:,i)./(P+1),wp(:,i)./(P+1),'color',ccolor_map(i,:),'linewidth',0.8....
                ,'LineStyle',char(line_style(i))....
                ,'Marker',char(marker_type(i))....
                ,'MarkerIndices',marker_indices,'MarkerSize',8....
                ),hold on
end
plot(K(:,1)./(P+1),K(:,1)./(P+1),':k')

xlabel('$K$','Interpreter','latex','FontSize',14);
ylabel('$\mathcal{R}$e $(K_{m})$'...
        ,'Interpreter','latex','FontSize',14)
xlim([0,pi])
xticks([0,pi/4,pi/2,3*pi/4,pi])
xticklabels({'\fontsize{14} 0','\fontsize{20}\pi\fontsize{14}/4'....
,'\fontsize{20}\pi\fontsize{14}/2'...
,'\fontsize{14}3\fontsize{20}\pi\fontsize{14}/4'...
,'\fontsize{20}\pi'})
ax.YAxis.FontSize=14;

legend({'$\beta=1.00$','$\beta=0.55$','$\beta=0.20$','$\beta=0.00$','exact'}...
    ,'FontSize',15,'Location','northwest','Interpreter','latex')

h.PaperUnits = 'inches';
h.PaperPosition = [0 0 6 4];
set(gca,'Color',[0.9 0.9 0.9]);

outerpos = ax.OuterPosition;
ti = ax.TightInset; 
%ti(4) = 0.01;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

if(plt_dump_flag==1)
    ffname = strcat('p',num2str(P),'_compare_beta_sdisc_wp');
    fname_png = strcat(png_fig,ffname,'.png');
    fname_eps = strcat(eps_fig,ffname,'.eps');
    print(fname_png,'-dpng','-r0')
    print(fname_eps,'-depsc','-r0')
end

% Dispersion error plot:
disper_err = zeros(length(K(:,1)),4);
h=figure;
ax = gca;
ax.TickLabelInterpreter = 'latex';
for i=1:4
    disper_err(:,i) = abs(wp(:,i) - K(:,i));
    marker_indices=[100:marker_step(i):length(K(:,i))];
    semilogy(K(:,i)./(P+1),disper_err(:,i)./(P+1),'color',ccolor_map(i,:)....
                ,'linewidth',0.8,'LineStyle',char(line_style(i))....
                ,'Marker',char(marker_type(i)).....
                ,'MarkerIndices',marker_indices,'MarkerSize',8),hold on
end

xx=xlabel('$K$','Interpreter','latex','FontSize',axes_labels_fontsize);
ylabel('$|\mathcal{R}$e$(K_{m})-K|$'....
    ,'Interpreter','latex','FontSize',axes_labels_fontsize)

xlim([0,pi/2])
ylim([1e-8,1])
grid on
xticks([0,pi/8,pi/4,3*pi/8,pi/2])
xticklabels({'\fontsize{14} 0'....
    '\fontsize{20}\pi\fontsize{14}/8'....
    ,'\fontsize{20}\pi\fontsize{14}/4'....
    ,'\fontsize{14}3\fontsize{20}\pi\fontsize{14}/8'....
    ,'\fontsize{20}\pi\fontsize{14}/2'})
ax.YAxis.FontSize=axes_ticks_fontsize;
legend(legend_wp,'FontSize',legend_fontsize,'Location'...
    ,'southeast','Interpreter','latex')

h.PaperUnits = 'inches';
h.PaperPosition = [0 0 7 5];
set(ax,'Color',[0.9 0.9 0.9]);

outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

if(plt_dump_flag==1)
    ffname = strcat('p',num2str(P),'_compare_beta_wp_err');
    fname_png = strcat(png_fig,ffname,'.png');
    fname_eps = strcat(eps_fig,ffname,'.eps');
    print(fname_png,'-dpng','-r0')
    print(fname_eps,'-depsc','-r0')
end

%------------------------------------
% Dissipation error plot:
%------------------------------------
dissip_err = zeros(length(K(:,1)),4);
h=figure;
ax = gca;
ax.TickLabelInterpreter = 'latex';
for i=1:4
    dissip_err(:,i) = abs(wd(:,i));
    marker_indices=[100:marker_step(i):length(K(:,i))];
    semilogy(K(:,i)./(P+1),dissip_err(:,i)./(P+1),'color',ccolor_map(i,:)....
                ,'linewidth',0.8,'LineStyle',char(line_style(i))....
                ,'Marker',char(marker_type(i)).....
                ,'MarkerIndices',marker_indices,'MarkerSize',8),hold on
end

xlabel('$K$','Interpreter','latex','FontSize',axes_labels_fontsize);
ylabel('$|\mathcal{I}$m$(K_{m})|$'....
    ,'Interpreter','latex','FontSize',axes_labels_fontsize)
xlim([0,pi/2])
ylim([1e-8,1])
grid on
xticks([0,pi/8,pi/4,3*pi/8,pi/2])
xticklabels({'\fontsize{14} 0'....
    '\fontsize{20}\pi\fontsize{14}/8'....
    ,'\fontsize{20}\pi\fontsize{14}/4'....
    ,'\fontsize{14}3\fontsize{20}\pi\fontsize{14}/8'....
    ,'\fontsize{20}\pi\fontsize{14}/2'})
ax.YAxis.FontSize=axes_ticks_fontsize;
legend(legend_wd,'FontSize',legend_fontsize,'Location'...
    ,'southeast','Interpreter','latex')

h.PaperUnits = 'inches';
h.PaperPosition = [0 0 7 5];
set(ax,'Color',[0.9 0.9 0.9]);

outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

if(plt_dump_flag==1)
    ffname = strcat('p',num2str(P),'_compare_beta_wd_err');
    fname_png = strcat(png_fig,ffname,'.png');
    fname_eps = strcat(eps_fig,ffname,'.eps');
    print(fname_png,'-dpng','-r0')
    print(fname_eps,'-depsc','-r0')
end

%========================
% total error plot:
%========================
error = zeros(length(K(:,1)),4);

h=figure;
ax = gca;
ax.TickLabelInterpreter = 'latex';
for i=1:4
    error(:,i) = abs(wd(:,i) - 1i.* wp(:,i) + 1i.*K(:,i));
    marker_indices=[100:marker_step(i):length(K(:,i))];
    semilogy(K(:,i)./(P+1),error(:,i)./(P+1),'color',ccolor_map(i,:)....
                ,'linewidth',0.8,'LineStyle',char(line_style(i))....
                ,'Marker',char(marker_type(i)).....
                ,'MarkerIndices',marker_indices,'MarkerSize',8),hold on
end

xlabel('$K$','Interpreter','latex','FontSize',axes_labels_fontsize);
ylabel('Total error'....
    ,'Interpreter','latex','FontSize',axes_labels_fontsize)

xlim([0,pi])
ylim([1e-6,10])
grid on
xticks([0,pi/4,pi/2,3*pi/4,pi])
xticklabels({'\fontsize{14} 0'....
    '\fontsize{20}\pi\fontsize{14}/4'....
    ,'\fontsize{20}\pi\fontsize{14}/2'....
    ,'\fontsize{14}3\fontsize{20}\pi\fontsize{14}/4'....
    ,'\fontsize{20}\pi'})
ax.YAxis.FontSize=axes_ticks_fontsize;
legend(legend_wp,'FontSize',legend_fontsize,'Location'...
    ,'southeast','Interpreter','latex')

h.PaperUnits = 'inches';
h.PaperPosition = [0 0 7 5];
set(ax,'Color',[0.9 0.9 0.9]);

outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4) -0.02;
ax.Position = [left bottom ax_width ax_height];

if(plt_dump_flag==1)
    ffname = strcat('p',num2str(P),'_compare_beta_total_err');
    fname_png = strcat(png_fig,ffname,'.png');
    fname_eps = strcat(eps_fig,ffname,'.eps');
    print(fname_png,'-dpng','-r0')
    print(fname_eps,'-depsc','-r0')
end