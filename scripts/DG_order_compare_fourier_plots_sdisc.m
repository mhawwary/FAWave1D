% Comparison of different DGp2 Beta Schemes:
close all
clear
clc

plt_dump_flag=1;

P=1:5;           % Polynomial order
Prk = 3;       % RK order

% Loading data:
Beta = 1.00;      % fully upwind flux
cfl_max=0.205;
CFL = 0.001*cfl_max;  % CFL = CFLmax


ccolor_map=hsv(length(P));  % Setting a colormap for Beta

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

% ccolor_map = 'black', 'blue' , 'cyan' , 'green' , 'red' 
ccolor_map = [0,0,0; 0,0,1; 0,1,1; 0,1,0; 1,0,1; 1,0,0;  0.5 1.0 0.5]; 
line_style={'-',':','--','-.','-','-','--',':'};
marker_type= {'o','^','*','s','<','+','v'};

% Dissipation plot:

h=figure;
ax = gca;
ax.TickLabelInterpreter = 'latex';

for i=1:length(P)
    clear K DG_wd DG_wp data
    
    main_dir = strcat(outdir,'DGp',num2str(P(i)),'/');
    dir = strcat(main_dir,'data/');

    fname_wd = strcat(dir,'DGp',num2str(P(i)).....
        ,'_Beta',num2str(Beta),'_sdisc_wd.dat');
    fname_wp = strcat(dir,'DGp',num2str(P(i))...
        ,'_Beta',num2str(Beta),'_sdisc_wp.dat');

    data=load(fname_wd);
   
    K = data(:,1)./(P(i)+1);
    if(P(i)==5)
        DG_wd = zeros(1,length(data(:,1)));
        ii = find(K>=2.066);
        jj = ii(1)-1;
        DG_wd(1:jj) = data(1:jj,6);
        DG_wd(jj+1:end) = data(jj+1:end,3);
        DG_wd(1,2)=0.5*(DG_wd(1,1)+DG_wd(1,3));
    elseif(P(i)==3)
        DG_wd= data(:,3);
        DG_wd(1,1)=DG_wd(2,1);
        DG_wd(end,1)=DG_wd(end-1,1);
    elseif(P(i)==4)
        DG_wd= data(:,5);
    else
        DG_wd= data(:,2);
        DG_wd(1,1)=DG_wd(2,1);
    end
    
    if(P(i)==1)
        marker_indices = [650,850,1050,length(K)];
    elseif(P(i)==2)
        marker_indices = [1100,1400,1650,length(K)];
    elseif(P(i)==3)
        marker_indices = [1450,1800,2250,length(K)];
    elseif(P(i)==4)
        length(K)
        marker_indices = [2050,2500,2850,length(K)];
    end    

    if(P(i)==5)
        plot(K,DG_wd./(P(i)+1),'color',ccolor_map(i,:),'linewidth',0.75....
            ,'LineStyle',char(line_style(i)))
    else
        marker_step = round(length(K)/5,0);
        plot(K,DG_wd./(P(i)+1),'color',ccolor_map(i,:),'linewidth',0.8....
            ,'LineStyle',char(line_style(i))....
            ,'Marker',char(marker_type(i))....
            ,'MarkerIndices',marker_indices,'MarkerSize',8....
            ,'MarkerFaceColor',ccolor_map(i,:)),hold on
    end
    
    disp(strcat('P= ',num2str(P(i))))

end
K_ex = linspace(0,pi,length(K));
plot(K_ex,0.0.*K_ex,':k')

xlabel('$K$','Interpreter','latex','FontSize',14);
ylabel('$\mathcal{I}$m$(K_{m})$'...
        ,'Interpreter','latex','FontSize',14)
% ylim([-40,10])
%yticks([-40:5:5])
xlim([0,pi])
xticks([0,pi/4,pi/2,3*pi/4,pi])
xticklabels({'\fontsize{14} 0','\fontsize{20}\pi\fontsize{14}/4'....
,'\fontsize{20}\pi\fontsize{14}/2'...
,'\fontsize{14}3\fontsize{20}\pi\fontsize{14}/4'...
,'\fontsize{20}\pi'})
ax.YAxis.FontSize=14;

legend({'P1','P2','P3','P4','P5','exact'}...
    ,'FontSize',15,'Location','southwest','Interpreter','latex')

h.PaperUnits = 'inches';
h.PaperPosition = [0 0 6 4];
set(gca,'Color',[0.9 0.9 0.9]);

outerpos = ax.OuterPosition;
ti = ax.TightInset; 
%ti(4) = 0.04;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

if(plt_dump_flag==1)
    ffname = strcat('DG_compareP_sdisc_wd');
    fname_png = strcat(outdir,png_fig,ffname);
    fname_eps = strcat(outdir,eps_fig,ffname);
    print(fname_png,'-dpng','-r0')
    print(fname_eps,'-depsc','-r0')
end

% Dispersion plot:

h=figure;
ax = gca;
ax.TickLabelInterpreter = 'latex';

for i=1:length(P)
    clear K DG_wd DG_wp data
    
    main_dir = strcat(outdir,'DGp',num2str(P(i)),'/');
    dir = strcat(main_dir,'data/');

    fname_wp = strcat(dir,'DGp',num2str(P(i))...
        ,'_Beta',num2str(Beta),'_sdisc_wp.dat');

    data=load(fname_wp);
   
    K = data(:,1)./(P(i)+1);
    if(P(i)==5)
        DG_wd = zeros(1,length(data(:,1)));
        ii = find(K>=2.066);
        jj = ii(1)-1;
        DG_wd(1:jj) = data(1:jj,6);
        DG_wd(jj+1:end) = data(jj+1:end,3);
        DG_wd(1,1)=0;
        DG_wd(1,2) = 0.5*(DG_wd(1,3)+DG_wd(1,1));
    elseif(P(i)==3)
        DG_wd= data(:,3);
        DG_wd(1,1)=DG_wd(2,1);
        DG_wd(end,1)=DG_wd(end-1,1);
    elseif(P(i)==4)
        DG_wd= data(:,5);
        DG_wd(1,1)=0;
        DG_wd(2,1) = 0.5*(DG_wd(3,1)+DG_wd(1,1));
    else
        DG_wd= data(:,2);
        DG_wd(1,1)=DG_wd(2,1);
    end

    if(P(i)==1)
        marker_indices = [750,950,1050,1150,length(K)];
    elseif(P(i)==2)
        marker_indices = [1000,1450,1650,length(K)];
    elseif(P(i)==3)
        marker_indices = [1425,2025,2250,length(K)];
    elseif(P(i)==4)
        length(K)
        marker_indices = [2000,2550,2800,length(K)];
    end    

    if(P(i)==5)
        plot(K,DG_wd./(P(i)+1),'color',ccolor_map(i,:),'linewidth',0.75....
            ,'LineStyle',char(line_style(i)))
    else
        marker_step = round(length(K)/5,0);
        plot(K,DG_wd./(P(i)+1),'color',ccolor_map(i,:),'linewidth',0.8....
            ,'LineStyle',char(line_style(i))....
            ,'Marker',char(marker_type(i))....
            ,'MarkerIndices',marker_indices,'MarkerSize',8....
            ,'MarkerFaceColor',ccolor_map(i,:)),hold on
    end
    
%     plot(K,DG_wd./(P(i)+1),'color',ccolor_map(i,:)....
%         ,'linewidth',1.5),hold on
    
    disp(strcat('P= ',num2str(P(i))))

end

plot(K,K,':k')

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

legend({'P1','P2','P3','P4','P5','exact'}...
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
    ffname = strcat('DG_compareP_sdisc_wp');
    fname_png = strcat(png_fig,ffname);
    fname_eps = strcat(eps_fig,ffname);
    print(fname_png,'-dpng','-r0')
    print(fname_eps,'-depsc','-r0')
end
% 
% % Dispersion error plot:
% disper_err(:,1) = abs(DG_wp1(:,1) - K(:,1));
% disper_err(:,2) = abs(DG_wp1(:,2) - K(:,2));
% disper_err(:,3) = abs(DG_wp1(:,3) - K(:,3));
% disper_err(:,4) = abs(DG_wp1(:,4) - K(:,4));
% 
% h=figure;
% ax = gca;
% ax.TickLabelInterpreter = 'latex';
% semilogy(K(:,1)./(P+1),disper_err(:,1),'-r','linewidth',1.5),hold on
% semilogy(K(:,2)./(P+1),disper_err(:,2),'-b','linewidth',1.5),hold on
% semilogy(K(:,3)./(P+1),disper_err(:,3),'-g','linewidth',1.5),hold on
% semilogy(K(:,4)./(P+1),disper_err(:,4),'-c','linewidth',1.5),hold on
% 
% xx=xlabel('K/P+1','Interpreter','latex','FontSize',15);
% yy=ylabel('Dispersion error','Interpreter','latex','FontSize',15);
% xlim([0,pi])
% xticks([0,pi/4,pi/2,3*pi/4,pi])
% xticklabels({'\fontsize{14} 0','\fontsize{20}\pi\fontsize{14}/4'....
% ,'\fontsize{20}\pi\fontsize{14}/2'...
% ,'\fontsize{14}3\fontsize{20}\pi\fontsize{14}/4'...
% ,'\fontsize{20}\pi'})
% ylim([1e-10,20])
% ax.YAxis.FontSize=14;
% set(ax, 'Units', 'Normalized');
% pos = get(ax, 'Position');
% offset = 0.03;
% set(ax,'Position',pos + [offset, offset, 0, 0])
% offset = 0.04;
% set(xx, 'Units', 'Normalized');
% pos = get(xx, 'Position');
% set(xx, 'Position', pos + [0, -offset, 0]);
% offset = 0.02;
% set(yy, 'Units', 'Normalized');
% pos = get(yy, 'Position');
% set(yy, 'Position', pos + [-offset, 0, 0]);
% 
% legend({'\beta=1.00','\beta=0.55','\beta=0.20','\beta=0.0'}....
%     ,'FontSize',16,'Location','southeast')
% set(gca,'Color',[0.9 0.9 0.9]);
% h.PaperUnits = 'inches';
% h.PaperPosition = [0 0 6 5];
% 
% fname = strcat(main_dir,'figures/p',num2str(P),'_compare_beta_sdisc_wp_error');
% print(fname,'-depsc','-r0')
% print(fname,'-dpng','-r0')
% 
% % Dissipation error plot:
% dissip_err(:,1) = abs(DG_wd1(:,1) - 0);
% dissip_err(:,2) = abs(DG_wd1(:,2) - 0);
% dissip_err(:,3) = abs(DG_wd1(:,3) - 0);
% 
% h=figure;
% ax = gca;
% ax.TickLabelInterpreter = 'latex';
% semilogy(K(:,1)./(P+1),dissip_err(:,1),'-r','linewidth',1.5),hold on
% semilogy(K(:,2)./(P+1),dissip_err(:,2),'-b','linewidth',1.5),hold on
% semilogy(K(:,3)./(P+1),dissip_err(:,3),'-g','linewidth',1.5),hold on
% 
% xx=xlabel('K/P+1','Interpreter','latex','FontSize',15);
% yy=ylabel('Dissipation error','Interpreter','latex','FontSize',15);
% xlim([0,pi])
% xticks([0,pi/4,pi/2,3*pi/4,pi])
% xticklabels({'\fontsize{14} 0','\fontsize{20}\pi\fontsize{14}/4'....
% ,'\fontsize{20}\pi\fontsize{14}/2'...
% ,'\fontsize{14}3\fontsize{20}\pi\fontsize{14}/4'...
% ,'\fontsize{20}\pi'})
% ylim([1e-10,20])
% ax.YAxis.FontSize=14;
% set(ax, 'Units', 'Normalized');
% pos = get(ax, 'Position');
% offset = 0.03;
% set(ax,'Position',pos + [offset, offset, 0, 0])
% offset = 0.04;
% set(xx, 'Units', 'Normalized');
% pos = get(xx, 'Position');
% set(xx, 'Position', pos + [0, -offset, 0]);
% offset = 0.02;
% set(yy, 'Units', 'Normalized');
% pos = get(yy, 'Position');
% set(yy, 'Position', pos + [-offset, 0, 0]);
% 
% legend({'\beta=1.00','\beta=0.55','\beta=0.20'}....
%     ,'FontSize',16,'Location','southeast')
% set(gca,'Color',[0.9 0.9 0.9]);
% h.PaperUnits = 'inches';
% h.PaperPosition = [0 0 6 5];
% 
% fname = strcat(main_dir,'figures/p',num2str(P),'_compare_beta_sdisc_wd_error');
% print(fname,'-depsc','-r0')
% print(fname,'-dpng','-r0')
% 
% % total error plot:
% error(:,1) = abs(DG_wd1(:,1) - 1i.* DG_wp1(:,1) + 1i.*K(:,1));
% error(:,2) = abs(DG_wd1(:,2) - 1i.* DG_wp1(:,2) + 1i.*K(:,1));
% error(:,3) = abs(DG_wd1(:,3) - 1i.* DG_wp1(:,3) + 1i.*K(:,1));
% error(:,4) = abs(DG_wd1(:,4) - 1i.* DG_wp1(:,4) + 1i.*K(:,1));
% 
% h=figure;
% ax = gca;
% ax.TickLabelInterpreter = 'latex';
% semilogy(K(:,1)./(P+1),error(:,1),'-r','linewidth',1.5),hold on
% semilogy(K(:,2)./(P+1),error(:,2),'-b','linewidth',1.5),hold on
% semilogy(K(:,3)./(P+1),error(:,3),'-g','linewidth',1.5),hold on
% semilogy(K(:,4)./(P+1),error(:,4),'-c','linewidth',1.5)
% 
% xx=xlabel('K/P+1','Interpreter','latex','FontSize',15);
% yy=ylabel('total error','Interpreter','latex','FontSize',15);
% xlim([0,pi])
% xticks([0,pi/4,pi/2,3*pi/4,pi])
% xticklabels({'\fontsize{14} 0','\fontsize{20}\pi\fontsize{14}/4'....
% ,'\fontsize{20}\pi\fontsize{14}/2'...
% ,'\fontsize{14}3\fontsize{20}\pi\fontsize{14}/4'...
% ,'\fontsize{20}\pi'})
% ylim([1e-10,20])
% ax.YAxis.FontSize=14;
% set(ax, 'Units', 'Normalized');
% pos = get(ax, 'Position');
% offset = 0.03;
% set(ax,'Position',pos + [offset, offset, 0, 0])
% offset = 0.04;
% set(xx, 'Units', 'Normalized');
% pos = get(xx, 'Position');
% set(xx, 'Position', pos + [0, -offset, 0]);
% offset = 0.02;
% set(yy, 'Units', 'Normalized');
% pos = get(yy, 'Position');
% set(yy, 'Position', pos + [-offset, 0, 0]);
% 
% legend({'\beta=1.00','\beta=0.55','\beta=0.20','\beta=0.00'}....
%     ,'FontSize',16,'Location','southeast')
% set(gca,'Color',[0.9 0.9 0.9]);
% h.PaperUnits = 'inches';
% h.PaperPosition = [0 0 6 5];
% 
% fname = strcat(main_dir,'figures/p',num2str(P),'_compare_beta_sdisc_error');
% print(fname,'-depsc','-r0')
% print(fname,'-dpng','-r0')
