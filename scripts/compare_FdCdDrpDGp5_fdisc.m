%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourier Stability plots for DGp5 vs CD6 vs FD6 and RK4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
close all
clc

plt_dump_flag=0;
ratio = 0.9;
nIter = 1.0;
approach = 'dt_const' ; % CFL_const

if(strcmp(approach,'dt_const'))
    sim_type = 'dtmax' ; % 'CFLmax'
elseif(strcmp(approach,'CFL_const'))
    sim_type = 'CFLmax' ; % 'CFLmax'
end

P = 5;      % Polynomial order
Prk = 4;         % RK order
cfl_max_dgu= 0.073;  % 0.0735
cfl_max_dgc = 0.103;  % 0.1035
cfl_max_dgb = 0.108; % beta =0.2
CFL_dgu = ratio * cfl_max_dgu;   % CFL no.
if(strcmp(approach,'CFL_const'))
    CFL_dgc = ratio * cfl_max_dgc ; 
    CFL_dgb = ratio * cfl_max_dgb ; 
elseif(strcmp(approach,'dt_const'))
    CFL_dgc = CFL_dgu;
    CFL_dgb = CFL_dgu;
end
% Beta = 0.0;      % upwind_parameter

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

legend_fontsize=17;
axes_labels_fontsize=16;
axes_ticks_fontsize=16;
marker_step=30;

%--------------------------------------------------------------------------
%                     DGp5 Loading or Computing data
%--------------------------------------------------------------------------

% load dominant mode data:
ss=strcat('loading DG upwind.........'); disp(ss);
[K_dgu,DGu_wd,DGu_wp]=load_dominant_mode(P,Prk,1,ratio); % upwind
ss=strcat('loading DG central.........'); disp(ss);
[K_dg,DGc_wd,DGc_wp]=load_dominant_mode(P,Prk,0,ratio); % central
ss=strcat('loading DG beta=',num2str(0.2),'.........'); disp(ss);
[K_dgb,DGb_wd,DGb_wp]=load_dominant_mode(P,Prk,0.2,ratio); % beta=0.2, RK4

G_dgc = exp(DGc_wd.*CFL_dgc.*nIter);
err_dgc_wd = abs(1-G_dgc);
err_dgc_wp = abs(K_dg - DGc_wp)./(P+1);

G_dgu = exp(DGu_wd.*CFL_dgu.*nIter);
err_dgu_wd = abs(1-G_dgu);
err_dgu_wp = abs(K_dgu - DGu_wp)./(P+1);
err_dgu_wd = abs(DGu_wd./(P+1));

G_dgb = exp(DGb_wd.*CFL_dgb.*nIter);
err_dgb_wd = abs(1-G_dgb);
err_dgb_wp = abs(K_dgb - DGb_wp)./(P+1);

kk = pi/4;
index(1) = find(K_dgu./(P+1) <= kk , 1, 'last' );

%--------------------------------------------------------------------------
% Explicit FD 6th order , central:
%--------------------------------------------------------------------------
ss=strcat('loading FD6 central.........'); disp(ss);
OA =6;
cfl_max_fd=[1.783,1.199];  % 1: is for cental
if(strcmp(approach, 'CFL_const'))
    CFL_fd(1) = ratio * cfl_max_fd(1);
    parent_dir = strcat('../../FD_FourierAnalysis_toolbox/');
    indir = strcat(parent_dir,'results/FD'...
        ,num2str(OA),'_RK',num2str(Prk),'/');
    % loadind data:
    fname= strcat(indir,datadir,'O',num2str(OA),'_RK',num2str(Prk)....
        ,'_',num2str(ratio),'CFLmax_wdwp.dat');
    data = load(fname);
    K_fd = data(:,1);
    wd_fd = data(:,2);
    wp_fd = data(:,3);
elseif(strcmp(approach, 'dt_const'))
    CFL_fd(1) = CFL_dgu.*(P+1);
    dK=0.005;
    K_fd = [0.005,dK:dK:pi-dK,pi];
    [~,FD]= FD_FourStab(OA,Prk, K_fd, 0,CFL_fd(1));
    K_fd= K_fd';
    wd_fd = FD.wd';
    wp_fd = FD.wp';
    clear FD
end

G_fd = exp(wd_fd.*CFL_fd(1).*nIter);
err_fd_wd = abs(wd_fd);
err_fd_wp = abs(K_fd - wp_fd);

index(2) = find(K_fd <= kk , 1, 'last' );

%--------------------------------------------------------------------------
% Explicit FD 6th order , 2-point upwind-biased :
%--------------------------------------------------------------------------
ss=strcat('loading FD6 upwind.........'); disp(ss);
OA =6;

if(strcmp(approach, 'CFL_const'))
    CFL_fd(2) = ratio * cfl_max_fd(2);
elseif(strcmp(approach, 'dt_const'))
    CFL_fd(2) = CFL_dgu.*(P+1);
end

[~,FD]= FD_FourStab(OA,Prk, K_fd, 2,CFL_fd(2));
wd_fd(:,2) = FD.wd;
wp_fd(:,2)= FD.wp;
clear FD

G_fd(:,2) = exp(wd_fd(:,2).*CFL_fd(2).*nIter);
err_fd_wd(:,2) = abs(1-G_fd(:,2));
err_fd_wd(:,2) = abs(wd_fd(:,2));
err_fd_wp(:,2) = abs(K_fd - wp_fd(:,2));

% index(2) = find(K_fd <= kk , 1, 'last' );

%--------------------------------------------------------------------------
% Compact CD 6th order :
%--------------------------------------------------------------------------
ss=strcat('loading CD6.........'); disp(ss);
OA =6;
CD_name= cell(1,3);
CD_name(1,1) = {strcat('CD$',num2str(OA))};

cfl_max_cd=1.421;
if(strcmp(approach, 'CFL_const'))
    CFL_cd = ratio * cfl_max_cd;
    parent_dir = strcat('../../FD_FourierAnalysis_toolbox/');
    indir = strcat(parent_dir,'results/CD'...
        ,num2str(OA),'_RK',num2str(Prk),'/');
    % loadind data:
    fname= strcat(indir,datadir,'O',num2str(OA),'_RK',num2str(Prk)....
        ,'_',num2str(ratio),'CFLmax_wdwp.dat');
    data = load(fname);
    K_cd = data(:,1);
    wd_cd = zeros(length(K_cd),3);
    G_cd = wd_cd; error_cd_wd = wd_cd;
    wd_cd(:,1) = data(:,2);
    wp_cd = data(:,3);
elseif(strcmp(approach, 'dt_const'))
    CFL_cd  = CFL_dgu.*(P+1);
    K_cd=K_fd;
    [~,CD]= CD_FourStab(OA,Prk, K_cd, CFL_cd);
    wd_cd = zeros(length(K_cd),3);
    G_cd = wd_cd; error_cd_wd = wd_cd; 
    wd_cd(:,1) = CD.wd(1,:)';
    wp_cd = CD.wp(1,:)';
    clear CD
end

G_cd(:,1) = exp(wd_cd(:,1).*CFL_cd.*nIter);
err_cd_wd(:,1) = abs(1-G_cd(:,1));
err_cd_wp = abs(K_cd - wp_cd);
%--------------------------
% Using PadeeFilters:
%--------------------------
alphaf = 0.49;
Of = 8;
CD_name(1,2) = {strcat('C$',num2str(OA)....
    ,'$F$',num2str(Of),'^{',num2str(alphaf),'}$')}; 
[G_cd(:,2)]=PadeeFilter(Of,alphaf,K_cd,G_cd(:,1));
wd_cd(1:end-1,2) = log(G_cd(1:end-1,2))./(CFL_cd.*nIter);
wd_cd(end,2)= wd_cd(end-1,2)-5;
err_cd_wd(:,2) = abs(1-G_cd(:,2));

alphaf = 0.40;
Of = 8;
CD_name(1,3) = {strcat('C$',num2str(OA)....
    ,'$F$',num2str(Of),'^{0.40}$')}; 
[G_cd(:,3)]=PadeeFilter(Of,alphaf,K_cd,G_cd(:,1));
wd_cd(1:end-1,3) = log(G_cd(1:end-1,3))./(CFL_cd.*nIter);
wd_cd(end,3)= wd_cd(end-1,3)-5;
err_cd_wd(:,3) = abs(1-G_cd(:,3));

err_cd_wd = abs(wd_cd);


index(3) = find(K_cd <= kk , 1, 'last' );
index1= index +1;


%--------------------------------------------------------------------------
% Explicit FDo11p of Bogey and Bailly:
%--------------------------------------------------------------------------
OA =2;
n_stages = 6;
RK_type = 'BogeyBailly6s';
RK_type = 'CLASSICAL';
Npts_FD = 7;
OA_f=6;
Npts_SF = 11; 
sigma_filter=1.0;
cfl_max_drp=2.053740;  % FDo11pSFo11p_RKo6s  
cfl_max_drp=1.891;     % FDo13pSFo13p_RKo6s
cfl_max_drp =  2.3980; %FDo7pSFo11p_RKo6s, DRP(4,7)
cfl_max_drp =  1.72; %FDo7pSFo11p_RK4
cfl_max_drp =  1.691; %Remez(2,7)_RK4
scheme_name = strcat('FDo',num2str(Npts_FD),'pSFo'....
    ,num2str(Npts_SF),'p_RKo',num2str(n_stages),'s');
drp_scheme_name = strcat('DRP($4$,$',num2str(Npts_FD),'$)-SF$'.....
    ,num2str(Npts_SF),'$-RK6s');
if(strcmp(RK_type,'CLASSICAL'))
    if(OA==2)
        drp_scheme_name = strcat('Remez($',num2str(OA),'$,$',num2str(Npts_FD),'$)-SF$'.....
        ,num2str(Npts_SF),'$');
    scheme_name = strcat('Remez',num2str(Npts_FD),'pSFo'....
        ,num2str(Npts_SF),'p_RK',num2str(Prk));
    else
    drp_scheme_name = strcat('DRP($',num2str(OA),'$,$',num2str(Npts_FD),'$)-SF$'.....
        ,num2str(Npts_SF),'$');
    scheme_name = strcat('FDo',num2str(Npts_FD),'pSFo'....
        ,num2str(Npts_SF),'p_RK',num2str(Prk));
    end

end
ss=strcat('constructing...',scheme_name,'.........'); disp(ss);
%scheme_name=strcat('DRP(4,',num2str(Npts_FD),')-SF(6,',num2str(Npts_SF),')-');

if(strcmp(approach, 'CFL_const'))
    CFL_drp = ratio * cfl_max_drp;
    parent_dir = strcat('../../FD_FourierAnalysis_toolbox/');
    indir = strcat(parent_dir,'results/',scheme_name,'/');
    % loadind data:
    fname= strcat(indir,datadir,scheme_name,'_',num2str(ratio)....
        ,'CFLmax_wdwp.mat');
    data = load(fname);
    K_drp = data.K;
    wd_drp = data.wd;
    wp_drp = data.wp;
elseif(strcmp(approach, 'dt_const'))
    CFL_drp = CFL_dgu.*(P+1);
    %CFL_drp = 0.05;
    dK=0.005;
    K_drp = [0.005,dK:dK:pi-dK,pi];
    [~,FD]= DRP_FourStab(OA,Npts_FD,RK_type,Prk,K_drp,CFL_drp);
    G_drp = exp(FD.wd.*CFL_drp.*nIter);
    G_drpsf=BogeyBaillyFilter(OA_f,Npts_SF,sigma_filter,K_drp,G_drp);
    wdtemp = log(G_drpsf)./CFL_drp;
    wd_drp = wdtemp';
    K_drp= K_drp';
    wp_drp = FD.wp';
    %wd_drp = FD.wd';
    clear FD
end

err_drp_wd = abs(wd_drp);
err_drp_wp = abs(K_drp - wp_drp);

%======================== dissip analysis calc ===========================%

K_test = [K_dgu(index(1))/(P+1), K_fd(index(2)), K_cd(index(3))];
K_test1 = [K_dgu(index1(1))/(P+1), K_fd(index1(2)), K_cd(index1(3))];

% DGp5:
slope = (G_dgu(index1(1)+1) - G_dgu(index(1)+1))/(K_test1(1)-K_test(1));
% ddK = kk-K_test(1);
% G_req(1) = slope * ddK +  G_dgu(index(1));
ddK = K_test1(1)-kk;
G_req(1) = G_dgu(index1(1)) - slope * ddK ;
% G_req(1) = G_dgu(index(1));

Ne_dg = 4;
Nv = [1,10];
dx_dg = 1.0/Ne_dg;
dx_fd = 1/((P+1)*Ne_dg);
a_dg = DGu_wp(index(1)) / ((P+1)*pi/4.0);
t_dg = Nv./a_dg
n_iter(1,:) = t_dg ./ (CFL_dgu * dx_dg);

% central FD:
a_fd = wp_fd(index(2),1) / (pi/4);
t_fdc = Nv./a_fd
n_iter(2,:) = t_fdc ./ (CFL_fd(1) * dx_fd);

% central FD:
a_fd = wp_fd(index(2),2) / (pi/4);
t_fd = Nv./a_fd
n_iter(3,:) = t_fd ./ (CFL_fd(2) * dx_fd);

a_cd = wp_cd(index(3)) / (pi/4);
t_cd = Nv./a_cd
n_iter(4,:) = t_cd ./ (CFL_cd * dx_fd);

round(n_iter)

% FD:
if(strcmp(approach,'CFL_const') && ratio==0.9)
    slope = (G_fd(index1(2),1) - G_fd(index(2),1))/(K_test1(2)-K_test(2));
    ddK = kk-K_test(2);
    G_req(2) = slope * ddK +  G_fd(index(2),1);
    slope = (G_fd(index1(2),2) - G_fd(index(2),2))/(K_test1(2)-K_test(2));
    ddK = kk-K_test(2);
    G_req(3) = slope * ddK +  G_fd(index(2),2);
elseif(strcmp(approach,'dt_const') && ratio==0.9)
    G_req(2) = G_fd(index(2),1);
    G_req(3) = G_fd(index(2),2);
end

% C6F8a0.4:
if(strcmp(approach,'CFL_const') && ratio==0.9)
    slope = (G_cd(index1(3),3) - G_cd(index(3),3))/(K_test1(3)-K_test(3));
    ddK = kk-K_test(3);
    G_req(4) = slope * ddK +  G_cd(index(3),3);
elseif(strcmp(approach,'dt_const') && ratio==0.9)
    G_req(4) = G_cd(index(3),3);
end

% C6F8a0.49:
if(strcmp(approach,'CFL_const') && ratio==0.9)
    slope = (G_cd(index1(3),2) - G_cd(index(3),2))/(K_test1(3)-K_test(3));
    ddK = kk-K_test(3);
    G_req(5) = slope * ddK +  G_cd(index(3),2);
elseif(strcmp(approach,'dt_const') && ratio==0.9)
    G_req(5) = G_cd(index(3),2);
end

if(strcmp(approach,'CFL_const') && ratio==0.9)
    err(1,1) = abs(1-G_req(1).^61);  err(1,2) = abs(1-G_req(1).^610); % DG
    err(2,1) = abs(1-G_req(2).^15);  err(2,2) = abs(1-G_req(2).^152); % FDc
    err(3,1) = abs(1-G_req(3).^22);  err(3,2) = abs(1-G_req(3).^224); % FDu
    err(4,1) = abs(1-G_req(4).^19);  err(4,2) = abs(1-G_req(4).^190); % a0.40
    err(5,1) = abs(1-G_req(5).^19);  err(5,2) = abs(1-G_req(5).^190); % a0.49
elseif(strcmp(approach,'dt_const') && ratio==0.9)
    err(1,1) = abs(1-G_req(1)^61);  err(1,2) = abs(1-G_req(1)^610); % DG
    err(2,1) = abs(1-G_req(2)^61);  err(2,2) = abs(1-G_req(2)^610); % FDc
    err(3,1) = abs(1-G_req(3)^61);  err(3,2) = abs(1-G_req(3)^608); % FDu
    err(4,1) = abs(1-G_req(4)^61);  err(4,2) = abs(1-G_req(4)^609); % a0.40
    err(5,1) = abs(1-G_req(5)^61);  err(5,2) = abs(1-G_req(5)^609); % a0.49
end



%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%                         Fully Discrete Plots                            %
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ss=strcat('G plot0.........'); disp(ss);

%-----------------------------------------------
% G plot0:
%-----------------------------------------------
h=figure;
ax = gca;
ax.TickLabelInterpreter = 'latex';
plot(K_dgu./(P+1),DGu_wd./(P+1),'-.k','linewidth',0.9),hold on
plot(K_cd,wd_cd(:,2),'->b'....
    ,'MarkerIndices',marker_step*0.5:marker_step:length(err_cd_wd(:,2))....
    ,'MarkerSize',8,'linewidth',0.9), hold on
plot(K_cd,wd_cd(:,3),'-sb'....
    ,'MarkerIndices',1:marker_step:length(err_cd_wd(:,3))....
    ,'MarkerSize',8,'linewidth',0.9), hold on
plot(K_fd,wd_fd(:,1),'--m'.....
    ,'MarkerIndices',1:marker_step:length(err_fd_wd(:,2))....
    ,'MarkerSize',8,'linewidth',0.9), hold on
plot(K_fd,wd_fd(:,2),'--oc'....
    ,'MarkerIndices',1:marker_step:length(err_fd_wd(:,2))....
    ,'MarkerSize',8,'linewidth',0.9), hold on
plot(K_drp,wd_drp,'--vg'....
    ,'MarkerIndices',1:marker_step:length(err_drp_wd)....
    ,'MarkerSize',8,'linewidth',0.9), hold on

plot(K_fd,0.0.*K_fd,':k')

xlabel('$K$','Interpreter','latex','FontSize',axes_labels_fontsize);
ylabel('$\mathcal{I}$m$(K_{m})$'...
    ,'Interpreter','latex','FontSize',axes_labels_fontsize)
xlim([0,pi])
ylim([-3.5,0.5])
grid on
xticks([0,pi/4,pi/2,3*pi/4,pi])
xticklabels({'\fontsize{14} 0'....
    ,'\fontsize{20}\pi\fontsize{14}/4'....
    ,'\fontsize{20}\pi\fontsize{14}/2'...
    ,'\fontsize{14}3\fontsize{20}\pi\fontsize{14}/4'...
    ,'\fontsize{20}\pi'})
ax.YAxis.FontSize=axes_ticks_fontsize;
legend_nn = {'DGp$5$-$\beta1.0$',char(CD_name(1,2))....
    ,char(CD_name(1,3)),'FD$6$-central','FD$6$-upwind',char(drp_scheme_name),'exact'};
legend(legend_nn,'FontSize',legend_fontsize,'Location','southwest','Interpreter','latex')
legend_nn = {'DGp$5$-$\beta1.0$',char(CD_name(1,2))....
    ,char(CD_name(1,3)),'FD$6$-central','FD$6$-upwind',char(drp_scheme_name)};

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
    ffname = strcat('p',num2str(P),'cd',num2str(OA),'_RK',num2str(Prk).....
        ,'_0',num2str(ratio*10),sim_type,'_G');
    fname_png = strcat(outdir,png_fig,ffname);
    fname_eps = strcat(outdir,eps_fig,ffname);
    print(fname_png,'-dpng','-r0')
    print(fname_eps,'-depsc','-r0')
end

%-----------------------------------------------
% Dissipation error plot:
%-----------------------------------------------
ss=strcat('dissip error plot.........'); disp(ss);
h=figure;
ax = gca;
ax.TickLabelInterpreter = 'latex';
% semilogy(K_dg./(P+1),err_dgc_wd,'-c','linewidth',1.5),hold on
% semilogy(K_dgb./(P+1),err_dgb_wd ,'-g','linewidth',1.5),hold on
semilogy(K_dgu./(P+1),err_dgu_wd ,'-.k','linewidth',0.9),hold on
semilogy(K_cd,err_cd_wd(:,2),'->b'....
    ,'MarkerIndices',marker_step*0.5:marker_step:length(err_cd_wd(:,2))....
    ,'MarkerSize',8,'linewidth',0.9),hold on
semilogy(K_cd,err_cd_wd(:,3),'-sb'....
    ,'MarkerIndices',1:marker_step:length(err_cd_wd(:,3))....
    ,'MarkerSize',9,'linewidth',0.9), hold on
semilogy(K_fd,err_fd_wd(:,1),'--m'.....
    ,'MarkerIndices',1:15:length(err_fd_wd(:,2))....
    ,'MarkerSize',8,'linewidth',0.9), hold on
semilogy(K_fd,err_fd_wd(:,2),'--oc'....
    ,'MarkerIndices',1:marker_step:length(err_fd_wd(:,2))....
    ,'MarkerSize',8,'linewidth',0.9), hold on
semilogy(K_drp,err_drp_wd,'--vg'....
    ,'MarkerIndices',1:marker_step:length(err_drp_wd)....
    ,'MarkerSize',8,'linewidth',0.9), hold on

xlabel('$K$','Interpreter','latex','FontSize',axes_labels_fontsize);
ylabel('$|\mathcal{I}$m$(K_{m})|$'....
    ,'Interpreter','latex','FontSize',axes_labels_fontsize)
xlim([0,pi/2])
ylim([1e-10,1])
grid on
xticks([0,pi/8,pi/4,3*pi/8,pi/2])
xticklabels({'\fontsize{14} 0'....
    '\fontsize{20}\pi\fontsize{14}/8'....
    ,'\fontsize{20}\pi\fontsize{14}/4'....
    ,'\fontsize{14}3\fontsize{20}\pi\fontsize{14}/8'....
    ,'\fontsize{20}\pi\fontsize{14}/2'})
ax.YAxis.FontSize=axes_ticks_fontsize;
legend(legend_nn,'FontSize',legend_fontsize....
    ,'Location','southeast','Interpreter','latex')

h.PaperUnits = 'inches';
h.PaperPosition = [0 0 8 6];
set(ax,'Color',[0.9 0.9 0.9]);

outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

if(plt_dump_flag==1)
    ffname = strcat('p',num2str(P),'cd',num2str(OA),'_RK',num2str(Prk).....
        ,'_0',num2str(ratio*10),sim_type,'_err_wd');
    fname_png = strcat(outdir,png_fig,ffname);
    fname_eps = strcat(outdir,eps_fig,ffname);
    print(fname_png,'-dpng','-r0')
    print(fname_eps,'-depsc','-r0')
end

%-----------------------------------------------
% Dispersion plot:
%-----------------------------------------------
ss=strcat('dispersion plot.........'); disp(ss);
h=figure;
ax = gca;
ax.TickLabelInterpreter = 'latex';
% plot(K_dg./(P+1),DGc_wp./(P+1),'-c','linewidth',1.5),hold on
% plot(K_dgb./(P+1),DGb_wp./(P+1),'-g','linewidth',1.5),hold on
plot(K_dgu./(P+1),DGu_wp./(P+1),'-.k','linewidth',0.9),hold on
plot(K_cd,wp_cd,'-b','linewidth',0.9), hold on
plot(K_fd,wp_fd(:,1),'--m'.....
    ,'MarkerIndices',1:marker_step:length(wp_fd(:,2))....
    ,'MarkerSize',8,'linewidth',0.9), hold on
plot(K_fd,wp_fd(:,2),'--oc'....
    ,'MarkerIndices',1:marker_step:length(wp_fd(:,2))....
    ,'MarkerSize',8,'linewidth',0.9), hold on
plot(K_drp,wp_drp,'--vg'....
    ,'MarkerIndices',1:marker_step:length(wp_drp)....
    ,'MarkerSize',8,'linewidth',0.9), hold on
plot(K_cd,K_cd,':k')

xlabel('$K$','Interpreter','latex','FontSize',axes_labels_fontsize);
ylabel('$\mathcal{R}$e $(K_{m})$'...
    ,'Interpreter','latex','FontSize',axes_labels_fontsize)

xlim([0,pi])
%     ylim([-20,1])
grid on
xticks([0,pi/4,pi/2,3*pi/4,pi])
xticklabels({'\fontsize{14} 0'....
    ,'\fontsize{20}\pi\fontsize{14}/4'....
    ,'\fontsize{20}\pi\fontsize{14}/2'...
    ,'\fontsize{14}3\fontsize{20}\pi\fontsize{14}/4'...
    ,'\fontsize{20}\pi'})
ax.YAxis.FontSize=axes_ticks_fontsize;
% legend({'DGp5, \beta=0.00','DGp5, \beta=0.20'....
%     ,'DGp5, \beta=1.00','CD6','FD6','exact'}....
%     ,'FontSize',14,'Location','northwest')

drp_scheme_name2= strcat('DRP($4$,$',num2str(Npts_FD),'$)-RK6s');
if(strcmp(RK_type,'CLASSICAL'))
    if(OA==2)
        drp_scheme_name2= strcat('Remez($',num2str(OA),'$,$',num2str(Npts_FD),'$)');
    else
        drp_scheme_name2= strcat('DRP($4$,$',num2str(Npts_FD),'$)');
    end
end
legend_nn = {'DGp$5$-$\beta1.0$','CD$6$','FD$6$-central','FD$6$-upwind'....
    ,drp_scheme_name2,'exact'};
legend(legend_nn,'FontSize',legend_fontsize,'Location','northwest','Interpreter','latex')

% set(ax, 'Units', 'Normalized');
% pos = get(ax, 'Position');
% offset = 0.04;
% set(ax,'Position',pos + [0, offset, 0, 0])
% offset = 0.04;
% set(xx, 'Units', 'Normalized');
% pos = get(xx, 'Position');
% set(xx, 'Position', pos + [0, -offset, 0]);
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
    ffname = strcat('p',num2str(P),'cd',num2str(OA),'_RK',num2str(Prk).....
        ,'_0',num2str(ratio*10),sim_type,'_wp');
    fname_png = strcat(outdir,png_fig,ffname);
    fname_eps = strcat(outdir,eps_fig,ffname);
    print(fname_png,'-dpng','-r0')
    print(fname_eps,'-depsc','-r0')
end

%-----------------------------------------------
% Dispersion error plot:
%-----------------------------------------------
ss=strcat('dispersion error plot.........'); disp(ss);
h=figure;
ax = gca;
ax.TickLabelInterpreter = 'latex';
% semilogy(K_dg./(P+1),err_dgc_wp,'-c','linewidth',1.5),hold on
% semilogy(K_dgb./(P+1),err_dgb_wp ,'-g','linewidth',1.5),hold on
semilogy(K_dgu./(P+1),err_dgu_wp ,'-.k','linewidth',0.9),hold on
semilogy(K_cd,err_cd_wp,'-b','linewidth',0.9), hold on
semilogy(K_fd,err_fd_wp(:,1),'--m'....
    ,'MarkerIndices',1:marker_step:length(err_fd_wp(:,1))....
    ,'linewidth',0.9), hold on
semilogy(K_fd,err_fd_wp(:,2),'--oc'....
    ,'MarkerIndices',1:marker_step:length(err_fd_wp(:,2))....
    ,'MarkerSize',8,'linewidth',0.9),hold on
semilogy(K_drp,err_drp_wp,'--vg'....
    ,'MarkerIndices',1:marker_step:length(err_drp_wp)....
    ,'MarkerSize',8,'linewidth',0.9), hold on

xx=xlabel('$K$','Interpreter','latex','FontSize',axes_labels_fontsize);
ylabel('$|\mathcal{R}$e$(K_{m})-K|$'....
    ,'Interpreter','latex','FontSize',axes_labels_fontsize)
xlim([0,pi/2])
ylim([1e-10,1])
grid on
xticks([0,pi/8,pi/4,3*pi/8,pi/2])
xticklabels({'\fontsize{14} 0'....
    '\fontsize{20}\pi\fontsize{14}/8'....
    ,'\fontsize{20}\pi\fontsize{14}/4'....
    ,'\fontsize{14}3\fontsize{20}\pi\fontsize{14}/8'....
    ,'\fontsize{20}\pi\fontsize{14}/2'})
ax.YAxis.FontSize=axes_ticks_fontsize;
% legend({'DGp5, \beta=0.00','DGp5, \beta=0.20'....
%     ,'DGp5, \beta=1.00','CD6','FD6'},'FontSize',14,'Location','southeast')
legend_nn = {'DGp$5$-$\beta1.0$','CD$6$','FD$6$-central','FD$6$-upwind',....
    drp_scheme_name2};
legend(legend_nn,'FontSize',legend_fontsize....
    ,'Location','southeast','Interpreter','latex')

% set(ax, 'Units', 'Normalized');
% pos = get(ax, 'Position');
% offset = 0.04;
% set(ax,'Position',pos + [0, offset, 0, 0])
% offset = 0.04;
% set(xx, 'Units', 'Normalized');
% pos = get(xx, 'Position');
% set(xx, 'Position', pos + [0, -offset, 0]);
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 8 6];
set(ax,'Color',[0.9 0.9 0.9]);

outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

if(plt_dump_flag==1)
    ffname = strcat('p',num2str(P),'cd',num2str(OA),'_RK',num2str(Prk).....
        ,'_0',num2str(ratio*10),sim_type,'_err_wp');
    fname_png = strcat(outdir,png_fig,ffname);
    fname_eps = strcat(outdir,eps_fig,ffname);
    print(fname_png,'-dpng','-r0')
    print(fname_eps,'-depsc','-r0')
end

