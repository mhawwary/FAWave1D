clear 
close all
clc

plt_dump_flag=0;
dK=0.001;
K=[0.001,2*dK:dK:pi-dK,pi];
OA = [6,6];
bias=[0,2]; 
bias_name = ['fully-upwind  ';'fully-upwind  ';'1-point biased'...
    ;'2-point biased';'1-point biased';'2-point biased';'1-point biased'];

line_style={':','--','-.','-','--','-','-.','-'};
marker_type= {'o','>','s','*','^','d','+','v','.'};

ccolor_map = [0,0,1; 0,0,0; 0,1,1; 1,0,1; 1,0,0; 0,1,0;....
    1.0 0.7 0.2; 0.5 0.2 0.8]; 
%ccolor_map = [0,0,0;0,0,0; 0,0,0; 0,0,0; 0,0,0; 0,0,0; 0,0,0; 0,0,0; 0,0,0; 0,0,0; 0,0,0;]; 
%ccolor_map=jet(length(OA));


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

marker_step=[400,450,400,450,500];
%===================Loading and Preparing schemes Data====================%

%==================================disper=================================%
h=figure;
ax = gca;
ax.TickLabelInterpreter = 'latex';
legend_temp = cell(1,6);
wd = zeros(8,length(K));
wp = wd;

% FD schemes:
for i=1:length(OA)
    [FD,~]= FD_FourStab(OA(i),3, K, bias(i), 1);
    wd(i,:)=FD.wd;
    wp(i,:) = FD.wp;
    marker_indices = 800:marker_step(i):length(K);
    plot(K,wp(i,:),'color',ccolor_map(i,:)...
    ,'linewidth',1.0,'LineStyle',char(line_style(i))....
    ,'Marker',char(marker_type(i))....
    ,'MarkerFaceColor',ccolor_map(i,:)....
    ,'MarkerIndices',marker_indices,'MarkerSize',8); hold on

    if(i==1)
        legend_temp(1,i) = { strcat('FD',num2str(OA(i)),', central')};
    else
        legend_temp(1,i) = { strcat('FD',num2str(OA(i)),', upwind-biased')};
    end
end

% DRP schemes:
% [BB,~]= DRP_FourStab(4,9,'CLASSICAL',0,K,0);  % BB(4,9)
[DRP,~]= DRP_FourStab(4,7,'CLASSICAL',0,K,0);  % DRP(4,7)
[Rem,~]= DRP_FourStab(2,7,'CLASSICAL',0,K,0);  % Remez(2,7)
wp(3,:) = DRP.wp;
wp(4,:) = Rem.wp;
G_temp = ones(1,length(K));
[G_bogey]=BogeyBaillyFilter(6,11,1.0,K,1.0);
wd(3,:)= 1-G_bogey;
marker_indices = 800:marker_step(4):length(K);
plot(K,wp(3,:),'color',ccolor_map(3,:)...
            ,'linewidth',1.0,'LineStyle',char(line_style(3))....
            ,'Marker',char(marker_type(3))....
            ,'MarkerFaceColor',ccolor_map(3,:)....
            ,'MarkerIndices',marker_indices,'MarkerSize',8); hold on
%legend_temp(1,3) = { strcat('DRP(4,7)')};
legend_temp(1,3) = { strcat('DRP$4$s$7$')};
marker_indices = 800:marker_step(5):length(K);
plot(K,wp(4,:),'color',ccolor_map(4,:)...
            ,'linewidth',1.0,'LineStyle',char(line_style(4))....
            ,'Marker',char(marker_type(4))....
            ,'MarkerFaceColor',ccolor_map(4,:)....
            ,'MarkerIndices',marker_indices,'MarkerSize',8); hold on
% legend_temp(1,4) = { strcat('Remez(2,7)')};
legend_temp(1,4) = { strcat('Rem$2$s$7$')};


% Compact schemes:
i=4; % alpha=0.4 filter
[G_cd]=PadeeFilter(8,0.40,K,1);
wd(i,:) = 1- G_cd;
i=5; % alpha = 0.49 filter
[G_cd]=PadeeFilter(8,0.49,K,1);
wd(i,:) = 1- G_cd;
%spatial scheme
i=5;
marker_indices = 800:marker_step(i-1):length(K);
[CD,~]= CD_FourStab(OA(1),3, K, 1);
wd(i,:)=CD.wd;
wp(i,:) = CD.wp;
plot(K,wp(i,:),'color',ccolor_map(i+1,:)...
        ,'linewidth',1.0,'LineStyle',char(line_style(i))....
        ,'Marker',char(marker_type(i))....
        ,'MarkerFaceColor',ccolor_map(i+1,:)....
        ,'MarkerIndices',marker_indices,'MarkerSize',8); hold on
legend_temp(1,i) = { strcat('CD',num2str(OA(1)))};

%DG scheme:
P=5;
beta=1;
ss=strcat('loading DGp',num2str(P),'-upwind.........'); disp(ss);
[K_dgu,wd_dgu,wp_dgu]= load_dominant_mode_mat(P,0,beta,0);
marker_indices = 100:30:length(K_dgu);
plot(K_dgu./(P+1),wp_dgu./(P+1),'color',ccolor_map(8,:)...
    ,'linewidth',1.0,'LineStyle',char(line_style(6))....
    ,'Marker',char(marker_type(6))....
    ,'MarkerFaceColor',ccolor_map(8,:)....
    ,'MarkerIndices',marker_indices,'MarkerSize',8), hold on
legend_temp(1,6) = {'DGp$5$-$\beta1.0$'};
plot(K,K,':k')


xx=xlabel('$K$','Interpreter','latex','FontSize',14);
yy=ylabel('$\mathcal{R}$e$(K_{m})$'....
    ,'Interpreter','latex','FontSize',14);
grid on
xlim([0,pi])
% ylim([-0.05,1.1])
xticks([0,pi/4,pi/2,3*pi/4,pi])
xticklabels({'\fontsize{14} 0','\fontsize{20}\pi\fontsize{14}/4'....
,'\fontsize{20}\pi\fontsize{14}/2'...
,'\fontsize{14}3\fontsize{20}\pi\fontsize{14}/4'...
,'\fontsize{20}\pi'})

legend(legend_temp,'FontSize',15....
    ,'Location','northwest','Interpreter','latex')

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
    ffname = strcat('DgFdCdDrp_compare_wp');
    fname_png = strcat(outdir,png_fig,ffname,'.png');
    fname_eps = strcat(outdir,eps_fig,ffname,'.eps');
    print(fname_png,'-dpng','-r0')
    print(fname_eps,'-depsc','-r0')
end

%================================disper_err=================================%
marker_step=[250,280,250,280,300];
h=figure;
ax = gca;
ax.TickLabelInterpreter = 'latex';
% FD schemes
for i=1:2
    error = abs(K-wp(i,:));
    marker_indices = 50:marker_step(i):length(K);
    semilogy(K,error,'color',ccolor_map(i,:)...
    ,'linewidth',1.0,'LineStyle',char(line_style(i))....
    ,'Marker',char(marker_type(i))....
    ,'MarkerFaceColor',ccolor_map(i,:)....
    ,'MarkerIndices',marker_indices,'MarkerSize',8); hold on
end

%DRP schemes
error = abs(K-wp(3,:));
marker_indices = 100:marker_step(1):length(K);
semilogy(K,error,'color',ccolor_map(3,:)...
            ,'linewidth',1.0,'LineStyle',char(line_style(3))....
            ,'Marker',char(marker_type(3)).....
            ,'MarkerFaceColor',ccolor_map(3,:)....
            ,'MarkerIndices',marker_indices,'MarkerSize',8); hold on
error = abs(K-wp(4,:));
marker_indices = 100:marker_step(5):length(K);
semilogy(K,error,'color',ccolor_map(4,:)...
            ,'linewidth',1.0,'LineStyle',char(line_style(4))....
            ,'Marker',char(marker_type(4))....
            ,'MarkerFaceColor',ccolor_map(4,:)....
            ,'MarkerIndices',marker_indices,'MarkerSize',8); hold on
% CD schemes
i=5;
error = abs(K-wp(i,:));
semilogy(K,error,'color',ccolor_map(i+1,:)...
        ,'linewidth',1.0,'LineStyle',char(line_style(i))....
        ,'Marker',char(marker_type(i))....
        ,'MarkerFaceColor',ccolor_map(i+1,:)....
        ,'MarkerIndices',marker_indices,'MarkerSize',8); hold on
% DG scheme:
error_dg = abs(K_dgu(:,1)-wp_dgu(:,1))./(P+1);
marker_indices = 5:15:length(K_dgu);
semilogy(K_dgu./(P+1),error_dg,'color',ccolor_map(8,:)...
    ,'linewidth',1.0,'LineStyle',char(line_style(6))....
    ,'Marker',char(marker_type(6))....
    ,'MarkerFaceColor',ccolor_map(8,:)....
    ,'MarkerIndices',marker_indices,'MarkerSize',8)

xx=xlabel('$K$','Interpreter','latex','FontSize',14);
yy=ylabel('$|\mathcal{R}$e$(K_{m})-K|$'....
    ,'Interpreter','latex','FontSize',14);
grid on
xlim([0,pi/2])
ylim([1e-8,1])
% xticks([0,pi/4,pi/2,3*pi/4,pi])
% xticklabels({'\fontsize{14} 0','\fontsize{20}\pi\fontsize{14}/4'....
% ,'\fontsize{20}\pi\fontsize{14}/2'...
% ,'\fontsize{14}3\fontsize{20}\pi\fontsize{14}/4'...
% ,'\fontsize{20}\pi'})
xticks([0,pi/8,pi/4,3*pi/8,pi/2])
xticklabels({'\fontsize{14} 0','\fontsize{20}\pi\fontsize{14}/8'....
,'\fontsize{20}\pi\fontsize{14}/4'...
,'\fontsize{14}3\fontsize{20}\pi\fontsize{14}/8'...
,'\fontsize{20}\pi\fontsize{14}/2'})

legend(legend_temp,'FontSize',15....
    ,'Location','northwest','Interpreter','latex')

ax.YAxis.FontSize=14;
set(gca,'Color',[0.9 0.9 0.9]);
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 7 5];

outerpos = ax.OuterPosition;
ti = ax.TightInset; 
ti(4)=0.05;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

if(plt_dump_flag==1)
    ffname = strcat('DgFdCdDrp_compare_wp_err');
    fname_png = strcat(outdir,png_fig,ffname,'.png');
    fname_eps = strcat(outdir,eps_fig,ffname,'.eps');
    print(fname_png,'-dpng','-r0')
    print(fname_eps,'-depsc','-r0')
end

%================================dissip_err=================================%
marker_step=[250,280,250,280,300];
h=figure;
ax = gca;
ax.TickLabelInterpreter = 'latex';
% FD schemes
for i=2:2
    error = abs(wd(i,:));
    marker_indices = 50:marker_step(i):length(K);
    semilogy(K,error,'color',ccolor_map(i,:)...
    ,'linewidth',1.0,'LineStyle',char(line_style(i))....
    ,'Marker',char(marker_type(i))....
    ,'MarkerFaceColor',ccolor_map(i,:)....
    ,'MarkerIndices',marker_indices,'MarkerSize',8); hold on
end

%DRP schemes
error = abs(wd(3,:));
marker_indices = 100:marker_step(1):length(K);
semilogy(K,error,'color',ccolor_map(3,:)...
            ,'linewidth',1.0,'LineStyle',char(line_style(3))....
            ,'Marker',char(marker_type(3)).....
            ,'MarkerFaceColor',ccolor_map(3,:)....
            ,'MarkerIndices',marker_indices,'MarkerSize',8); hold on

% CD schemes
i=4;
error = abs(wd(i,:));
marker_indices = 100:marker_step(5):length(K);
semilogy(K,error,'color',ccolor_map(4,:)...
            ,'linewidth',1.0,'LineStyle',char(line_style(1))....
            ,'Marker',char(marker_type(1))....
            ,'MarkerFaceColor',ccolor_map(4,:)....
            ,'MarkerIndices',marker_indices,'MarkerSize',8); hold on
    
i=5;
error = abs(wd(i,:));
semilogy(K,error,'color',ccolor_map(i+1,:)...
        ,'linewidth',1.0,'LineStyle',char(line_style(i))....
        ,'Marker',char(marker_type(i))....
        ,'MarkerFaceColor',ccolor_map(i+1,:)....
        ,'MarkerIndices',marker_indices,'MarkerSize',8); hold on
% DG scheme:
error_dg = abs(wd_dgu(:,1))./(P+1);
marker_indices = 5:15:length(K_dgu);
semilogy(K_dgu./(P+1),error_dg,'color',ccolor_map(8,:)...
    ,'linewidth',1.0,'LineStyle',char(line_style(6))....
    ,'Marker',char(marker_type(6))....
    ,'MarkerFaceColor',ccolor_map(8,:)....
    ,'MarkerIndices',marker_indices,'MarkerSize',8)

xx=xlabel('$K$','Interpreter','latex','FontSize',14);
yy=ylabel('$|\mathcal{I}$m$(K_{m})|$'....
    ,'Interpreter','latex','FontSize',14);
grid on
xlim([0,pi/2])
ylim([1e-8,1])
% xticks([0,pi/4,pi/2,3*pi/4,pi])
% xticklabels({'\fontsize{14} 0','\fontsize{20}\pi\fontsize{14}/4'....
% ,'\fontsize{20}\pi\fontsize{14}/2'...
% ,'\fontsize{14}3\fontsize{20}\pi\fontsize{14}/4'...
% ,'\fontsize{20}\pi'})
xticks([0,pi/8,pi/4,3*pi/8,pi/2])
xticklabels({'\fontsize{14} 0','\fontsize{20}\pi\fontsize{14}/8'....
,'\fontsize{20}\pi\fontsize{14}/4'...
,'\fontsize{14}3\fontsize{20}\pi\fontsize{14}/8'...
,'\fontsize{20}\pi\fontsize{14}/2'})

legend({'FD$6$, upwind-biased', 'SF$11$','CF$^{0.40}$'....
    ,'CF$^{0.49}$','DGp$5$-$\beta1.0$'},'FontSize',15....
    ,'Location','northwest','Interpreter','latex')

ax.YAxis.FontSize=14;
set(gca,'Color',[0.9 0.9 0.9]);
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 7 5];

outerpos = ax.OuterPosition;
ti = ax.TightInset; 
ti(4)=0.05;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

if(plt_dump_flag==1)
    ffname = strcat('DgFdCdDrp_compare_wd_err');
    fname_png = strcat(outdir,png_fig,ffname,'.png');
    fname_eps = strcat(outdir,eps_fig,ffname,'.eps');
    print(fname_png,'-dpng','-r0')
    print(fname_eps,'-depsc','-r0')
end
