%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to dump dominant modes by calling dump_dominant_mode function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc

plt_dump_flag=0;
dump_flag=1;
dump_ii=2;
option_cent= 2;

P=3;           % Polynomial order
Prk = 3;       % RK order
Beta = 1.0;      % upwind_parameter
CFL_max = 0.130;

ratio = 0.4;
CFL = ratio * CFL_max;

%====================== Loaading the data ================================%

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

% figdir='figures/';
% datadir = 'data/';
% indir=strcat('./results/','DGp',num2str(P),'_RK',num2str(Prk),'/');
% fname = strcat(indir,datadir,'DGp',num2str(P),'_RK',num2str(Prk)....
%     ,'_Beta',num2str(Beta),'_',num2str(ratio),'CFLmax_');
% fname_wd = strcat(fname,'wd.dat');
% fname_wp = strcat(fname,'wp.dat');
% data_wd = load(fname_wd);
% data_wp = load(fname_wp);
% K = data_wd(:,1);
% Nk = length(K);
% wd = zeros(Nk,P+1);
% wp=wd; 
% 
% for i=1:P+1
%     wd(:,i) = data_wd(:,i+1);
%     wp(:,i) = data_wp(:,i+1);
% end


% wd(1:113,1) = wd(1:113,2);
% wd(190:end,1) = wd(190:end,3);
% wd(239:end,1) = wd(239:end,4);
% wp(1:113,1) = wp(1:113,2);
% wp(190:end,1) = wp(190:end,3);
% wp(239:end,1) = wp(239:end,4);

if(option_cent==1)
    %This is for option(1) of central flux and all other Beta's
    K = [linspace(0.01,(P+1)*pi-0.1,50*(P+1)),(P+1)*pi]; % for true solution
%     dK=0.2;
%     K=[0.001,0.002,dK:dK:(P+1)*pi-dK,(P+1)*pi];
    [~,DG]= DG_FourStab(P,Prk, K, Beta, CFL,1.0);
else
    %This is for option(2) of central flux and all other Beta's
    dK=0.01;
    Kp =[0.001,0.002,dK:dK:(P+1)*pi-dK,(P+1)*pi];
    %Kn = flip(-Kp);
    %K = [Kn,Kp];
    K=Kp;
    [~,DG]= DG_FourStab(P,Prk, K, Beta, CFL,3);
end

marker_type= {'o','^','s','*','<','+','v'};
marker_indices = 1:1:length(K);

Nk = length(K);
wd = zeros(Nk,P+1);
wp=wd; 

wd_name = cell(1,P+1);
wp_name = wd_name;

for i=1:P+1
    wd_name(1,i) = {strcat('wd', num2str(i))};
    wp_name(1,i) = {strcat('wp', num2str(i))};   
end

for i=1:P+1
   wd(:,i) = DG.(char(wd_name(1,i)))';
   wp(:,i) = DG.(char(wp_name(1,i)))';
end

% for p3-RK3-beta1, ratio=0.5:
temp = wp(end,4);
wp(end,4)=wp(end,2);
wp(end,2) = temp;
temp = wd(end,4);
wd(end,4)=wd(end,2);
wd(end,2) = temp;

% for p3-RK4-beta1:
% temp = wp(end,1);
% wp(end,1)=wp(end,2);
% wp(end,2) = temp;
% temp = wd(end,1);
% wd(end,1)=wd(end,2);
% wd(end,2) = temp;

% for p3-RK4-beta0.2:
% ktest = 1.405*(P+1);
% [wd,wp]=switch_phymode_curves(K,dK,wd,wp,2,1,ktest);
% ktest = 2.98*(P+1);
% [wd,wp]=switch_phymode_curves(K,dK,wd,wp,4,1,ktest);


%======================  Dissipation curves  =============================%
ccolor_map = [0,0,1; 0,1,1; 0,0,0; 1,0,1; 0,1,0; 1,0,0]; 
h=figure;
ax = gca;
ax.TickLabelInterpreter = 'latex';
for i=1:P+1
    plot(K./(P+1),wd(:,i),'color',ccolor_map(i,:)...
            ,'linewidth',1.5,'LineStyle','-'), hold on
end

xlim([0,pi])
% ylim([-0.16,0.06])
xticks([0,pi/4,pi/2,3*pi/4,pi])
xticklabels({'\fontsize{14} 0','\fontsize{20}\pi\fontsize{14}/4'....
,'\fontsize{20}\pi\fontsize{14}/2'...
,'\fontsize{14}3\fontsize{20}\pi\fontsize{14}/4'...
,'\fontsize{20}\pi'})
grid on
xx=xlabel('K/P+1','Interpreter','latex','FontSize',15);
yy=ylabel('$\mathcal{I}$m(K$_{m}$)','Interpreter','latex','FontSize',15);
legend({'mode$_{1}$','mode$_{2}$','mode$_{3}$','mode$_{4}$'},'FontSize',15....
    ,'Location','northwest','Interpreter','latex');
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
% 
if(plt_dump_flag==1)
    ffname = strcat('p',num2str(P),'_RK',num2str(Prk).....
        ,'_beta',num2str(Beta),'_0',num2str(ratio*10),'CFLmax_wd');
    fname_png = strcat(png_fig,ffname);
    fname_eps = strcat(eps_fig,ffname);
    print(fname_png,'-dpng','-r0')
    print(fname_eps,'-depsc','-r0')
end

%======================  Dispersion curves  ==============================%
h=figure;
ax = gca;
ax.TickLabelInterpreter = 'latex';
for i=1:P+1
    plot(K./(P+1),wp(:,i)./(P+1),'color',ccolor_map(i,:)...
            ,'linewidth',1.5,'LineStyle','-'....
            ,'Marker',char(marker_type(i))....
            ,'MarkerIndices',marker_indices,'MarkerSize',8), hold on
end
plot(K./(P+1),K./(P+1),':k','linewidth',1.0)

xlim([0,pi])
% ylim([-3,5])
xticks([0,pi/4,pi/2,3*pi/4,pi])
xticklabels({'\fontsize{14} 0','\fontsize{20}\pi\fontsize{14}/4'....
,'\fontsize{20}\pi\fontsize{14}/2'...
,'\fontsize{14}3\fontsize{20}\pi\fontsize{14}/4'...
,'\fontsize{20}\pi'})
grid on
xx=xlabel('K/P+1','Interpreter','latex','FontSize',15);
yy=ylabel('$\mathcal{R}$e(K$_{m}$)','Interpreter','latex','FontSize',15);
legend({'mode$_{1}$','mode$_{2}$','mode$_{3}$','mode$_{4}$','exact'},'FontSize',15....
    ,'Location','northeast','Interpreter','latex');
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
        ,'_beta',num2str(Beta),'_0',num2str(ratio*10),'CFLmax_wp');
    fname_png = strcat(png_fig,ffname);
    fname_eps = strcat(eps_fig,ffname);
    print(fname_png,'-dpng','-r0')
    print(fname_eps,'-depsc','-r0')
end

if(dump_flag==1)
    disp('dumping......')
    dump_dominant_mode_mat(P,Prk,Beta,ratio,K',wd(:,dump_ii),wp(:,dump_ii));
end



