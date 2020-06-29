clear
close all
clc
%--------------------------------------------------------------
% Script to plot DGp1  CFL for stability with different Beta's
%---------------------------------------------------------------
true_tol = 1.0;

P=5;           % Polynomial order
Prk = 4;       % RK order
%Beta = 0.0:0.05:1.0;      % upwind_parameter
Beta = [0.0,1.0];

CFL_test=0.1:0.001:0.105;
K_test = 0.005:0.005:(P+1)*pi;

CFL_max = zeros(1,length(Beta));
cfl_K = CFL_max;

for j=1:length(Beta)
    if (j==2 )
        CFL_test=0.07:0.001:0.075;
    end
    [CFL_max(j),cfl_K(j)]= identify_stable_CFL_DG(P,Prk, K_test, Beta(j), CFL_test,1e-10);
end
% 
% % Saving the output:
%Preparing directories:
%------------------------
resdir='../results/';
outdir=strcat(resdir,'DGp',num2str(P),'_RK',num2str(Prk),'/');
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

% %figdir = './matlab_results/';
% %dir = './data/';
% fname = strcat(outdir,'DGp',num2str(P),'_RK',num2str(Prk),'_stability_limits.dat');
% fileID = fopen(fname,'w');
% print_data= [Beta;CFL_max];
% fprintf(fileID,'%1.2f %1.3f\r\n',print_data);
% fclose(fileID);
% 
% 
% % Plotting:
% 
[max_cfl,ii] = max(CFL_max);
% 
% fprintf('Beta with Max CFL is, Beta:%1.2f, CFL:%1.3f\n',Beta(ii),max_cfl);
% 
% h=figure;
% 
% plot(Beta,CFL_max,'--or','linewidth',1.5)
% xlabel('\bf\beta','FontSize',15)
% ylabel('\bfCFL_{max}','FontSize',15)
% 
% dtitle = strcat('\bfStability limits for DGp'....
%     ,num2str(P),' and RK',num2str(Prk),' vs upwind parameter \beta');
% title(dtitle,'FontSize',14);
% 
% set(gca,'Color',[0.8 0.8 0.8]);
% 
% figname = strcat(outdir,'DGp',num2str(P),'_RK',num2str(Prk)....
%         ,'_stability_limits');
% saveas(h,figname,'png');
% savefig(h,figname,'compact');


% fname = strcat(outdir,'DGp',num2str(P),'_RK'....
%     ,num2str(Prk),'_stability_limits.dat');
% 
% data=load(fname);
% bbeta= data(:,1);
% cfl__max= data(:,2);
% 
% Prk=3;
% fname1 = strcat(outdir,'DGp',num2str(P),'_RK'....
%     ,num2str(Prk),'_stability_limits.dat');
% 
% clear data
% 
% data=load(fname1);
% bbeta1= data(:,1);
% cfl__max1= data(:,2);
% 
% Prk=4;
% fname2 = strcat(outdir,'DGp',num2str(P),'_RK'....
%     ,num2str(Prk),'_stability_limits.dat');
% 
% clear data
% 
% data=load(fname2);
% bbeta2= data(:,1);
% cfl__max2= data(:,2);
% 
% h=figure;
% 
% plot(bbeta,cfl__max,'-ok','linewidth',1.5....
%     ,'markerfacecolor','r','markeredgecolor','r','markersize',8),hold on
% plot(bbeta1,cfl__max1,'-ob','linewidth',1.5....
%     ,'markerfacecolor','g','markeredgecolor','g','markersize',8),hold on
% plot(bbeta2,cfl__max2,'--^r','linewidth',1.5....
%     ,'markerfacecolor','k','markeredgecolor','k','markersize',8),hold on
% 
% 
% xlabel('\bf\beta','FontSize',15)
% ylabel('\bfCFL_{max}','FontSize',15)
% 
% dtitle = strcat('\bfStability limits for DGp'....
%     ,num2str(P),' vs upwind parameter \beta');
% 
% title(dtitle,'FontSize',14);
% 
% legend({'RK2','RK3','RK4'},'fontsize',15)
% 
% set(gca,'Color',[0.8 0.8 0.8]);
% 
% figname = strcat(outdir,'DGp',num2str(P),'_comapreRK',num2str(Prk)....
%         ,'_stability_limits');
% saveas(h,figname,'png');
% savefig(h,figname,'compact');

