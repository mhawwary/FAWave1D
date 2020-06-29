%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to dump dominant modes by calling dump_dominant_mode function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc

plt_dump_flag=0;
dump_flag=1;
dump_ii=1;
option_cent= 1;

P=2;           % Polynomial order
Beta = 0.55;      % upwind_parameter
Prk=0;

if(option_cent==1)
    %This is for option(1) of central flux and all other Beta's
    K = [linspace(0.01,(P+1)*pi-0.1,50*(P+1)),(P+1)*pi]; % for true solution
    dK=0.1;
    K=[0.001,0.002,dK:dK:(P+1)*pi-dK,(P+1)*pi];
    [DG,~]= DG_FourStab(P,Prk, K, Beta, 0,2.0);
else
    %This is for option(2) of central flux and all other Beta's
    dK=0.20;
    Kp =[0.001,0.002,dK:dK:(P+1)*pi-dK,(P+1)*pi];
    Kn = flip(-Kp);
    K = [Kn,Kp];
    [DG,~]= DG_FourStab(P,Prk, K, Beta, 1,3.0);
end

marker_type= {'o','^','s','*','<','+','v'};
marker_indices = 20:1:length(K);
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


%for p3-beta1:
%temp = wd(end,3);
%wd(end,3)=wd(end,1);
%wd(end,1)=temp;

% for p3-beta0 & beta0.2:
% ktest = 1.4*(P+1);
% ix = find(abs(K-ktest)<dK);
% kerr = abs(K(ix)-ktest);
% [~,ii] = min(kerr);
% index = ix(ii);
% temp = wp(index:end,1);
% wp(index:end,1) = wp(index:end,2);
% wp(index:end,2) = temp;
% temp = wd(index:end,1);
% wd(index:end,1) = wd(index:end,2);
% wd(index:end,2) = temp;
% 
% ktest = 3*(P+1);
% ix = find(abs(K-ktest)<dK);
% kerr = abs(K(ix)-ktest);
% [~,ii] = min(kerr);
% index = ix(ii);
% temp = wp(index:end,2);
% wp(index:end,2) = wp(index:end,4);
% wp(index:end,4) = temp;
% temp = wd(index:end,2);
% wd(index:end,2) = wd(index:end,4);
% wd(index:end,4) = temp;


% wd(1:113,1) = wd(1:113,2);
% wd(190:end,1) = wd(190:end,3);
% wd(239:end,1) = wd(239:end,4);
% wp(1:113,1) = wp(1:113,2);
% wp(190:end,1) = wp(190:end,3);
% wp(239:end,1) = wp(239:end,4);

%--------------DGp5
% temp = wd(:,1);
% wd(:,1) = wd(:,3);
% wd(:,3) = temp;
% temp = wp(:,1);
% wp(:,1) = wp(:,3);
% wp(:,3) = temp;
% 
% temp = wd(110,1);
% wd(110,1) = wd(110,6);
% wd(110,6) = temp;
% temp = wp(110,1);
% wp(110,1) = wp(110,6);
% wp(110,6) = temp;
% temp = wd(111:117,1);
% wd(111:117,1) = wd(111:117,5);
% wd(111:117,5) = temp;
% temp = wp(111:117,1);
% wp(111:117,1) = wp(111:117,5);
% wp(111:117,5) = temp;
%------------------------
% 
% temp = wd(end,1);
% wd(end,1) = wd(end,2);
% wd(end,2)=temp;
% temp = wp(end,1);
% wp(end,1) = wp(end,2);
% wp(end,2)=temp;

%======================  Dissipation curves  =============================%
ccolor_map = [0,0,1; 0,1,1; 0,0,0; 1,0,1; 0,1,0; 1,0,0]; 
h=figure;
ax = gca;
ax.TickLabelInterpreter = 'latex';
for i=1:P+1
    plot(K./(P+1),wd(:,i),'color',ccolor_map(i,:)...
            ,'linewidth',1.5,'LineStyle','-'), hold on
end

xlim([-pi,pi])
% ylim([-3,5])
xticks([-pi,-3*pi/4,-pi/2,-pi/4,0,pi/4,pi/2,3*pi/4,pi])
xticklabels({'\fontsize{14}-\fontsize{20}\pi'....
    ,'\fontsize{14}-3\fontsize{20}\pi\fontsize{14}/4'....
    ,'\fontsize{14}-\fontsize{20}\pi\fontsize{14}/2'....
    ,'\fontsize{14}-\fontsize{20}\pi\fontsize{14}/4'....
    ,'\fontsize{14} 0','\fontsize{20}\pi\fontsize{14}/4'....
,'\fontsize{20}\pi\fontsize{14}/2'...
,'\fontsize{14}3\fontsize{20}\pi\fontsize{14}/4'...
,'\fontsize{20}\pi'})
grid on
xx=xlabel('$K$','Interpreter','latex','FontSize',15);
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

xlim([-pi,pi])
% ylim([-3,5])
xticks([-pi,-3*pi/4,-pi/2,-pi/4,0,pi/4,pi/2,3*pi/4,pi])
xticklabels({'\fontsize{14}-\fontsize{20}\pi'....
    ,'\fontsize{14}-3\fontsize{20}\pi\fontsize{14}/4'....
    ,'\fontsize{14}-\fontsize{20}\pi\fontsize{14}/2'....
    ,'\fontsize{14}-\fontsize{20}\pi\fontsize{14}/4'....
    ,'\fontsize{14} 0','\fontsize{20}\pi\fontsize{14}/4'....
,'\fontsize{20}\pi\fontsize{14}/2'...
,'\fontsize{14}3\fontsize{20}\pi\fontsize{14}/4'...
,'\fontsize{20}\pi'})
grid on
xx=xlabel('$K$','Interpreter','latex','FontSize',15);
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

%============================  Dumping =================================%

if(dump_flag==1)
    disp('dumping......')
    dump_dominant_mode_mat(P,0,Beta,0,K',wd(:,dump_ii),wp(:,dump_ii))
end


% clear
% close all
% clc
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % A script to dump and plot semi-discete dispersion/dissipation of DG
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_flag =1;
% plt_dump_flag=0;
% dump_flag=0;
% 
% true_tol = 1.0;
% P= 1:5;           % Polynomial order
% Prk = 0;       % RK order
% Beta = 1.00;      % upwind_parameter
% CFL=0.1;
% 
% ccolor_map=hot(length(P));  % Setting a colormap for Beta
% 
% dK=0.005;
% 
% figdir='figures/';
% datadir = 'data/';
% png_fig = 'figures/png/';
% eps_fig = 'figures/eps/';
% if(Beta==0)
%     beta_name = '0_00';
% elseif(Beta==1)
%     beta_name = '1_00';
% else
%     beta_name = strcat('0',num2str(100*Beta));
% end
% 
% wd_name = cell(1,6);
% wp_name = wd_name;
% 
% for i=1:6
%     
%     wd_name(1,i) = {strcat('wd', num2str(i))};
%     wp_name(1,i) = {strcat('wp', num2str(i))};
%     
% end
% 
% 
% if(dump_flag==1)
%     for i=1:length(P)
%         K=[0.00:dK:(P(i)+1)*pi-dK,(P(i)+1)*pi];
%         [DGsd,~]= DG_FourStab(P(i),Prk, K, Beta, CFL,true_tol);
%         
%         outdir=strcat('./results/DGp',num2str(P(i)),'/');
% 
%         if(dump_flag==1)
%             Nk = length(K);
%             print_data =zeros(Nk,P(i)+2);
%             print_data(:,1)=K;
%             ss= strcat('saving output.....p= ',num2str(P(i)));
%             disp(ss);
%             % Saving the output:
%             fname = strcat(outdir,datadir,'DGp',num2str(P(i)),....
%                 '_Beta',num2str(Beta),'_sdisc_wd.dat');
%             fileID = fopen(fname,'w');
% 
%             ff = '%1.3f';
%             for j=1:P(i)+1
%                 print_data(:,j+1)= DGsd.(char(wd_name(1,j)))(1,:);
%                 ff = strcat(ff,' %1.10e');
%             end
%             ff = strcat(ff,'\r\n');
%             for m=1:Nk
%                 fprintf(fileID,ff,print_data(m,:));
%             end
%             fclose(fileID);
%             %clear fname print_data
% 
%             fname = strcat(outdir,datadir,'DGp',num2str(P(i)),....
%                 '_Beta',num2str(Beta),'_sdisc_wp.dat');
%             fileID = fopen(fname,'w');
% 
%             ff = '%1.3f';
%             for j=1:P(i)+1
%                 print_data(:,j+1)= DGsd.(char(wp_name(1,j)))(1,:);
%                 ff = strcat(ff,' %1.10e');
%             end
%             ff = strcat(ff,'\r\n');
%             for m=1:Nk
%                 fprintf(fileID,ff,print_data(m,:));
%             end
%             fclose(fileID);
% 
%     %         clear K print_data
% 
%         end
% 
%         if(plot_flag==1)
%             if(P(i)==1)
%                 plot(K./(P(i)+1),DGsd.wd1(1,:)./(P(i)+1)), hold on
%             else
%                 plot(K./(P(i)+1),DGsd.wd3(1,:)./(P(i)+1)), hold on
%             end
% 
%         end
% 
%     end
% end




