clear
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A script to dump and plot semi-discete dispersion/dissipation of DG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_flag =1;
plt_dump_flag=0;
dump_flag=0;

true_tol = 1.0;
P= 2;           % Polynomial order
Prk = 0;       % RK order
Beta = 0.00;      % upwind_parameter
CFL=0.1;

ccolor_map=hot(length(P));  % Setting a colormap for Beta

dK=0.005;

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

if(Beta==0)
    beta_name = '0_00';
elseif(Beta==1)
    beta_name = '1_00';
else
    beta_name = strcat('0',num2str(100*Beta));
end

%==========================================================================
%                  Single order semi-discrete plots
%==========================================================================
if(plot_flag ==1)
    %This is for option(1) of central flux
    dK=0.0025;
    K=[-(P+1)*pi,-(P+1)*pi+dK:dK:-dK,-0.002,-0.001,0.001,0.002,dK:dK:(P+1)*pi-dK,(P+1)*pi];
    [DGsd,~]= DG_FourStab(P,0, K, Beta, 0,1.0);
    %This is for option(2) of central flux and all other Beta's
%     dK=0.20;
%     Kp =[0.001,0.002,dK:dK:(P+1)*pi-dK,(P+1)*pi];
%     Kn = flip(-Kp);
%     K = [Kn,Kp];
%     [DGsd,~]= DG_FourStab(P,Prk, K, Beta, 1,3.0);
    sum1 = abs(DGsd.we1(1,:)).^2 + abs(DGsd.we2(1,:)).^2....
        +abs(DGsd.we3(1,:)).^2;
    we1 = 100.*abs(DGsd.we1(1,:)).^2./sum1;
    we2 = 100.*abs(DGsd.we2(1,:)).^2./sum1;
    we3 = 100.*abs(DGsd.we3(1,:)).^2./sum1;
    
%     we3(7543:7545) = [100,100,100];
%     we5(1) = 0.0016;
%     we1(1) = 97.21;
    %outdir=strcat(outdir,'DGp',num2str(P),'/');
    marker_step = 500;
    marker_indices = 1:marker_step:length(K);
    %======================= Dissipation plot =====================%
    h=figure;
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    plot(K./(P+1),DGsd.wd1(1,:)./(P+1),'-sk','linewidth',0.9),hold on
    plot(K./(P+1),DGsd.wd2(1,:)./(P+1),'--vk','linewidth',0.9),hold on
    plot(K./(P+1),DGsd.wd3(1,:)./(P+1),':k','linewidth',0.9),hold on
    %     plot(K./(p+1),zeros(1,length(K)),'--k')

    xlabel('$K$','Interpreter','latex','FontSize',14);
    ylabel('$\mathcal{I}$m$(K_{m})$'...
        ,'Interpreter','latex','FontSize',14)
    ylim([-4,1.0])
    xlim([-pi,pi])
    xticks([-pi,-3*pi/4,-pi/2,-pi/4,0,pi/4,pi/2,3*pi/4,pi])
    xticklabels({'-\fontsize{20}\pi'....
        ,'-\fontsize{14}3\fontsize{20}\pi\fontsize{14}/4'....
        ,'-\fontsize{20}\pi\fontsize{14}/2'....
        ,'-\fontsize{20}\pi\fontsize{14}/4','\fontsize{14} 0'....
        ,'\fontsize{20}\pi\fontsize{14}/4'....
        ,'\fontsize{20}\pi\fontsize{14}/2'...
        ,'\fontsize{14}3\fontsize{20}\pi\fontsize{14}/4'...
        ,'\fontsize{20}\pi'})
    ax.YAxis.FontSize=14;
    legend({'mode$(1)$','mode$(2)$','mode$(3)$'}....
        ,'FontSize',15,'Location','northwest','Orientation','horizontal',....
        'Interpreter','latex')

    h.PaperUnits = 'inches';
    h.PaperPosition = [0 0 8 5.5];
    set(gca,'Color',[0.9 0.9 0.9]);

    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    ti(4) = 0.05;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];

    if(plt_dump_flag==1)
        ffname = strcat('p',num2str(P),'_beta'...
            ,beta_name,'_sdisc_wd');
        fname_png = strcat(png_fig,ffname,'.png');
        fname_eps = strcat(eps_fig,ffname,'.eps');
        print(fname_png,'-dpng','-r0')
        print(fname_eps,'-depsc','-r0')
    end

    %======================= Dispersion plot =====================%
    h=figure;
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    plot(K./(P+1),DGsd.wp3(1,:)./(P+1),'-sk','MarkerIndices',marker_indices,'MarkerSize',8,'linewidth',0.9),hold on
    plot(K./(P+1),DGsd.wp1(1,:)./(P+1),'--vk','MarkerIndices',marker_indices,'MarkerSize',8,'linewidth',0.9),hold on
    plot(K./(P+1),DGsd.wp2(1,:)./(P+1),'-.ok','MarkerIndices',marker_indices,'MarkerSize',8,'linewidth',0.9),hold on
    plot(K./(P+1),K./(P+1),':k','linewidth',0.70)

    xlabel('$K$','Interpreter','latex','FontSize',14);
    ylabel('$\mathcal{R}$e $(K_{m})$'...
        ,'Interpreter','latex','FontSize',14)
    h.PaperUnits = 'inches';
    h.PaperPosition = [0 0 8 5.5];
    set(gca,'Color',[0.9 0.9 0.9]);

    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    ti(4) = 0.05;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];

    xlim([-pi,pi])
    ylim([-4,4])
    xticks([-pi,-3*pi/4,-pi/2,-pi/4,0,pi/4,pi/2,3*pi/4,pi])
    xticklabels({'-\fontsize{20}\pi'....
        ,'-\fontsize{14}3\fontsize{20}\pi\fontsize{14}/4'....
        ,'-\fontsize{20}\pi\fontsize{14}/2'....
        ,'-\fontsize{20}\pi\fontsize{14}/4','\fontsize{14} 0'....
        ,'\fontsize{20}\pi\fontsize{14}/4'....
        ,'\fontsize{20}\pi\fontsize{14}/2'...
        ,'\fontsize{14}3\fontsize{20}\pi\fontsize{14}/4'...
        ,'\fontsize{20}\pi'})
    ax.YAxis.FontSize=14;
    legend({'mode$(1)$','mode$(2)$','mode$(3)$','exact'}....
        ,'FontSize',15,'Location','northwest','Orientation','horizontal',....
        'Interpreter','latex')
    
    if(plt_dump_flag==1)
        ffname = strcat('p',num2str(P),'_beta'...
            ,beta_name,'_sdisc_wp');
        fname_png = strcat(png_fig,ffname,'.png');
        fname_eps = strcat(eps_fig,ffname,'.eps');
        print(fname_png,'-dpng','-r0')
        print(fname_eps,'-depsc','-r0')
    end
    
    %%%%%%%%%%%%%%%%%%%% Weights plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    h=figure;
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    plot(K./(P+1), we3,'-k','linewidth',1.0),hold on
    plot(K./(P+1), we2,'--k','linewidth',1.0),hold on
    plot(K./(P+1), we1,':k','linewidth',1.0), hold on
    
    xlabel('K','Interpreter','latex','FontSize',15);
    ylabel('mode weight in $\% $'.....
        ,'Interpreter','latex','FontSize',15)
%     title(strcat('BR2p',num2str(p),', $\eta=$',num2str(eta)),'Interpreter'...
%         ,'latex','FontSize',15)
%     xlim([0,(P+1)*pi])
%     xticks([0:pi:(P+1)*pi])
%     ax.XTickLabel={'\fontsize{14} 0'....
%     ,'\fontsize{20}\pi'....
%     ,'\fontsize{14}2\fontsize{20}\pi'...
%     ,'\fontsize{14}3\fontsize{20}\pi'...
%     ,'\fontsize{14}4\fontsize{20}\pi'....
%     ,'\fontsize{14}5\fontsize{20}\pi'};
 xlim([0,pi])
    xticks([0,pi/4,pi/2,3*pi/4,pi])
    ax.XTickLabel={'\fontsize{14} 0'....
    ,'\fontsize{20}\pi/\fontsize{14}4'....
    ,'\fontsize{20}\pi/\fontsize{14}2'...
    ,'\fontsize{14}3\fontsize{20}\pi/\fontsize{14}4'...
    ,'\fontsize{20}\pi'};
    ylim([0,100])

    yyaxis right
    ax1 = gca;
    right_color = [0.5,0.8,1.0];
    right_color=[0.5,0.5,1.0];
    right_color1=[0,1,0];
    ax1.YColor = right_color;

    plot(K./(P+1) ,DGsd.wp3(1,:)./(P+1),'Marker','s','color',right_color,'LineStyle'...
        ,'-','MarkerIndices',marker_indices,'MarkerSize',8....
        ,'linewidth',0.5),hold on
    plot(K./(P+1) ,DGsd.wp2(1,:)./(P+1),'Marker','o','color',right_color,'LineStyle'....
        ,'--','MarkerIndices',marker_indices,'MarkerSize',8....
        ,'linewidth',0.5),hold on
    plot(K./(P+1) ,DGsd.wp1(1,:)./(P+1),'Marker','v','color',right_color,'LineStyle'....
        ,':','MarkerIndices',marker_indices,'MarkerSize',8....
        ,'linewidth',0.5),hold on
    % set(ax1, 'Box', 'off'); 
    %ylim([-10,10])
    ylabel('$\mathcal{R}$e ($K_{m}$)','Interpreter'...
        ,'latex','FontSize',15)
    % legend({'\alpha_{1}','\alpha_{2}','\lambda_{1}','\lambda_{2}'}....
    %     ,'FontSize',15,'Location','northeastoutside')
    legend({'mode(1)','mode(2)','mode(3)'}....
        ,'FontSize',15,'Location','northeastoutside')


    set(ax,'Color',[0.9 0.9 0.9]);
    h.PaperUnits = 'inches';
    h.PaperPosition = [0 0 10.5 5.5];

    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3)+0.2;
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];

    if(plt_dump_flag==2)
        ffname = strcat('p',num2str(P),'_beta'...
            ,beta_name,'_sdisc_weights');
        fname_png = strcat(png_fig,ffname,'.png');
        fname_eps = strcat(eps_fig,ffname,'.eps');
        print(fname_png,'-dpng','-r0')
        print(fname_eps,'-depsc','-r0')
    end

    
end
