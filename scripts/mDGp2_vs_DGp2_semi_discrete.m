
clear
close all
clc
%------------------------------------------
% Script to plot DGp2 fourier foots
%------------------------------------------
true_tol = 1.0;

P=2;           % Polynomial order
Prk = 3;       % RK order

% bb = 0.0:0.05:1.0;
% Nb = length(bb); 
% aa = 0.0:0.005:1.0;
% Na = length(aa);
N=P+1;

Beta = ones(P+1,1);      % upwind_parameter
%Alpha = ones(P+1,1);      % upwind_parameter
Alpha=[1;1.00;0.21];
%Beta(end,:) = bb;

CFL_max=0.62;
CFL=[0.005,(0.1:0.1:1.00)*CFL_max]; 
K = 0.005:0.005:(P+1)*pi;

ccolor_map=hot(length(CFL));  % Setting a colormap for Beta

if (Beta==1)
    flux = ' fully upwind flux';
elseif(Beta==0)
    flux = ' central flux';
else
    flux = strcat(' hybrid flux (\beta= ',num2str(Beta),')');
end

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%                         Fully Discrete Plots                            %
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%-------------------------------------------------------------------
%      fully disc Dissipation
%-------------------------------------------------------------------
h1=figure;

[mDGsd,~]= mDG_FourStab(P,Prk, K, Beta, Alpha, 1,true_tol);
[DGsd,~]= DG_FourStab(P,Prk, K, 1, 1,true_tol);
[DGsd_b,~]= DG_FourStab(P,Prk, K, 0.55, 1,true_tol);
plot(K, mDGsd.wd1(1,:),'-r','linewidth',1.5),hold on
plot(K, DGsd.wd1(1,:),'-b','linewidth',1.5),hold on
plot(K, DGsd_b.wd1(1,:),'-k','linewidth',1.5),hold on
  
xlabel('K (wave number)','FontSize',15)
ylabel('\Omega_{imag}','FontSize',15)
xlim([0,(P+1)*pi])
xticks(0:pi:(P+1)*pi)
xticklabels({'\fontsize{14} 0','\fontsize{22}\pi'....
    ,'\fontsize{14}2\fontsize{22}\pi'....
    ,'\fontsize{14}3\fontsize{22}\pi'})

legend({'fully upwind mDGp2, \alpha_{p}=0.21','fully upwind DGp2'....
    ,'hybrid DGp2, \beta=0.55'},'FontSize',15)

dtitle = strcat('Semi-Discrete Dissipation analysis, mDGp'....
    ,num2str(P),' vs DGp',num2str(P),' with  ',flux.....
    ,'. \alpha=[',num2str(Alpha(1)),', '....
    ,num2str(Alpha(2)),', ',num2str(Alpha(3)), '], jump parameter');
title(dtitle,'FontSize',14);

set(gca,'Color',[0.8 0.8 0.8]);


%-------------------------------------------------------------------
%      fully disc Dispersion
%-------------------------------------------------------------------
h=figure;

plot(K, mDGsd.wp1(1,:),'-r','linewidth',1.5),hold on
plot(K, DGsd.wp1(1,:),'-b','linewidth',1.5),hold on
plot(K, DGsd_b.wp1(1,:),'-k','linewidth',1.5),hold on
     
plot(K,K,'--k')

xlabel('K (wave number)','FontSize',15)
ylabel('\Omega_{real}','FontSize',15)
xlim([0,(P+1)*pi])
xticks(0:pi:(P+1)*pi)
xticklabels({'\fontsize{14} 0','\fontsize{22}\pi'....
    ,'\fontsize{14}2\fontsize{22}\pi'....
    ,'\fontsize{14}3\fontsize{22}\pi'})

legend({'fully upwind mDGp2, \alpha_{p}=0.21','fully upwind DGp2'....
    ,'hybrid DGp2, \beta=0.55'},'FontSize',15)

ptitle = strcat('Semi-Discrete Dispersion analysis, mDGp'....
    ,num2str(P),' vs DGp',num2str(P),' with  ',flux.....
    ,'. \alpha=[',num2str(Alpha(1)),', '....
    ,num2str(Alpha(2)),', ',num2str(Alpha(3)), '], jump parameter');
title(ptitle,'FontSize',14);

set(gca,'Color',[0.8 0.8 0.8]);

clear figname

