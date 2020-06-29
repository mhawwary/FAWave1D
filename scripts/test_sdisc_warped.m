clear
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A script to dump and plot semi-discete dispersion/dissipation of DG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_flag =1;
plt_dump_flag=1;
dump_flag=0;

true_tol = 1.0;
P= 4;           % Polynomial order
Prk = 0;       % RK order
Beta = 1.0;      % upwind_parameter
CFL=0.0;
X=-0.26;

ccolor_map=hot(length(P));  % Setting a colormap for Beta

dK=0.005;

figdir='figures/';
datadir = 'data/';
png_fig = 'figures/png/';
eps_fig = 'figures/eps/';
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
    K=[-(P+1)*pi,-(P+1)*pi+dK:dK:-dK,-0.002,-0.001,0.001,0.002.....
        ,dK:dK:(P+1)*pi-dK,(P+1)*pi];
    [DGsd,~]=  DG_FourStab_Warped(P,Prk, K, Beta,0.0,true_tol,X);
    %This is for option(2) of central flux and all other Beta's
%     dK=0.20;
%     Kp =[0.001,0.002,dK:dK:(P+1)*pi-dK,(P+1)*pi];
%     Kn = flip(-Kp);
%     K = [Kn,Kp];
%     [DGsd,~]= DG_FourStab(P,Prk, K, Beta, 1,3.0);
%     sum1 = abs(DGsd.we1(1,:)).^2 + abs(DGsd.we2(1,:)).^2....
%         +abs(DGsd.we3(1,:)).^2;
    sum1 = abs(DGsd.we1(1,:)).^2 + abs(DGsd.we2(1,:)).^2....
        +abs(DGsd.we3(1,:)).^2+abs(DGsd.we4(1,:)).^2....
        +abs(DGsd.we5(1,:)).^2;
    we1 = 100.*abs(DGsd.we1(1,:)).^2./sum1;
    we2 = 100.*abs(DGsd.we2(1,:)).^2./sum1;
    we3 = 100.*abs(DGsd.we3(1,:)).^2./sum1;
    we4 = 100.*abs(DGsd.we4(1,:)).^2./sum1;
    we5 = 100.*abs(DGsd.we5(1,:)).^2./sum1;
    
%     for m=1:length(K)
%        wee = DGsd.we1(m).*DGsd.V1(:,m); 
%        we1(m) = norm(wee); 
%        wee = DGsd.we2(m).*DGsd.V2(:,m); 
%        we2(m) = norm(wee);
%        wee = DGsd.we3(m).*DGsd.V3(:,m); 
%        we3(m) = norm(wee);
%        sum1(m) = abs(we1(1,m))^.2 + abs(we2(1,m)).^2....
%         +abs(we3(1,m)).^2;
%         we1(m) = 100.*abs(we1(1,m)).^2./sum1(m);
%         we2(m) = 100.*abs(we2(1,m)).^2./sum1(m);
%         we3(m) = 100.*abs(we3(1,m)).^2./sum1(m);
%     end
    
    outdir=strcat('./results/DGp',num2str(P),'/');
    marker_step = 550;
    marker_indices = 1:marker_step:length(K);
    %======================= Dissipation plot =====================%
    h=figure;
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    plot(K./(P+1),DGsd.wd1(1,:)./(P+1),'-sk'....
        ,'MarkerIndices',marker_indices,'MarkerSize',8....
        ,'linewidth',2),hold on
    plot(K./(P+1),DGsd.wd2(1,:)./(P+1),'--vb'....
        ,'MarkerIndices',marker_indices,'MarkerSize',8....
        ,'linewidth',2),hold on
    plot(K./(P+1),DGsd.wd3(1,:)./(P+1),':om'....
        ,'MarkerIndices',marker_indices,'MarkerSize',8....
        ,'linewidth',2),hold on
    plot(K./(P+1),DGsd.wd4(1,:)./(P+1),'-.+c'....
        ,'MarkerIndices',marker_indices,'MarkerSize',8....
        ,'linewidth',2),hold on
    plot(K./(P+1),DGsd.wd5(1,:)./(P+1),'-r'....
        ,'MarkerIndices',marker_indices,'MarkerSize',8....
        ,'linewidth',2),hold on
    %     plot(K./(p+1),zeros(1,length(K)),'--k')

    xlabel('$K$','Interpreter','latex','FontSize',14);
    ylabel('$\mathcal{I}$m$(K_{m})$'...
        ,'Interpreter','latex','FontSize',14)
%     ylim([-4,1.0])
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
    legend({'mode$(1)$','mode$(2)$','mode$(3)$,','mode$(4)$','mode$(5)$'}....
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
            ,beta_name,'_sdisc_wd_Xw',num2str(X));
        fname_png = strcat(outdir,png_fig,ffname,'.png');
        fname_eps = strcat(outdir,eps_fig,ffname,'.eps');
        print(fname_png,'-dpng','-r0')
        print(fname_eps,'-depsc','-r0')
    end

    %======================= Dispersion plot =====================%
    h=figure;
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    plot(K./(P+1),DGsd.wp1(1,:)./(P+1),'-sk'....
        ,'MarkerIndices',marker_indices,'MarkerSize',8....
        ,'linewidth',2),hold on
    plot(K./(P+1),DGsd.wp2(1,:)./(P+1),'--vb'....
        ,'MarkerIndices',marker_indices,'MarkerSize',8....
        ,'linewidth',2),hold on
    plot(K./(P+1),DGsd.wp3(1,:)./(P+1),':om'....
        ,'MarkerIndices',marker_indices,'MarkerSize',8....
        ,'linewidth',2),hold on
    plot(K./(P+1),DGsd.wp4(1,:)./(P+1),'-.+c'....
        ,'MarkerIndices',marker_indices,'MarkerSize',8....
        ,'linewidth',2),hold on
    plot(K./(P+1),DGsd.wp5(1,:)./(P+1),'-r'....
        ,'MarkerIndices',marker_indices,'MarkerSize',8....
        ,'linewidth',2),hold on
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
%     ylim([-4,4])
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
            ,beta_name,'_sdisc_wp_Xw',num2str(X));
        fname_png = strcat(outdir,png_fig,ffname,'.png');
        fname_eps = strcat(outdir,eps_fig,ffname,'.eps');
        print(fname_png,'-dpng','-r0')
        print(fname_eps,'-depsc','-r0')
    end
    
    %%%%%%%%%%%%%%%%%%%% Weights plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    h=figure;
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    plot(K./(P+1), we1,'-k','linewidth',2.0),hold on
    plot(K./(P+1), we2,'--b','linewidth',2.0),hold on
    plot(K./(P+1), we3,':m','linewidth',2.0), hold on
    plot(K./(P+1), we4,'-.c','linewidth',2.0),hold on
    plot(K./(P+1), we5,'-r','linewidth',2.0)

    xlabel('K','Interpreter','latex','FontSize',15);
    ylabel('mode weight in $\% $'.....
        ,'Interpreter','latex','FontSize',15)
%     title(strcat('BR2p',num2str(p),', $\eta=$',num2str(eta)),'Interpreter'...
%         ,'latex','FontSize',15)
    xlim([0,pi])
    xticks([0,pi/4,pi/2,3*pi/4,pi])
    xticklabels({'\fontsize{14} 0'....
        ,'\fontsize{20}\pi\fontsize{14}/4'....
        ,'\fontsize{20}\pi\fontsize{14}/2'...
        ,'\fontsize{14}3\fontsize{20}\pi\fontsize{14}/4'...
        ,'\fontsize{20}\pi'})
    ylim([0,100])

%     yyaxis right
%     ax1 = gca;
%     right_color = [0.5,0.8,1.0];
%     right_color=[0.5,0.5,1.0];
%     ax1.YColor = right_color;
% 
%     plot(K ,DGsd.wp1(1,:),'Marker','s','color',right_color,'LineStyle'...
%         ,'-','MarkerIndices',marker_indices,'MarkerSize',8....
%         ,'linewidth',0.5),hold on
%     plot(K ,DGsd.wp2(1,:),'Marker','v','color',right_color,'LineStyle'....
%         ,'--','MarkerIndices',marker_indices,'MarkerSize',8....
%         ,'linewidth',0.5),hold on
%     plot(K ,DGsd.wp3(1,:),'Marker','o','color',right_color,'LineStyle'....
%         ,':','MarkerIndices',marker_indices,'MarkerSize',8....
%         ,'linewidth',0.5),hold on
%   
%     % set(ax1, 'Box', 'off'); 
%     ylim([-10,10])
%     ylabel('$\mathcal{R}$e ($K_{m}$)','Interpreter'...
%         ,'latex','FontSize',15)
%     % legend({'\alpha_{1}','\alpha_{2}','\lambda_{1}','\lambda_{2}'}....
%     %     ,'FontSize',15,'Location','northeastoutside')
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

    if(plt_dump_flag==1)
        ffname = strcat('p',num2str(P),'_beta'...
            ,beta_name,'_sdisc_weights_Xw',num2str(X));
        fname_png = strcat(outdir,png_fig,ffname,'.png');
        fname_eps = strcat(outdir,eps_fig,ffname,'.eps');
        print(fname_png,'-dpng','-r0')
        print(fname_eps,'-depsc','-r0')
    end

    
end
