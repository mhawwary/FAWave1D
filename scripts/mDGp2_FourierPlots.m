
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

for j=1:length(CFL)
    
     [~,DGfd]= mDG_FourStab(P,Prk, K, Beta, Alpha, CFL(j),true_tol);
     plot(K, DGfd.wd1(1,:),'color',ccolor_map(j,:),'linewidth',1.5),hold on
     
    % Saving the output:
%     fname = strcat(dir,'DGp',num2str(P),'_RK',num2str(Prk)....
%         ,'_Beta',num2str(Beta),'_CFL',num2str(CFL(j)),'_wd.dat');
%     fileID = fopen(fname,'w');
%     print_data= [K;DGfd.wd1(1,:);DGfd.wd2(1,:)];
%     fprintf(fileID,'%1.3f %1.10e %1.10e\r\n',print_data);
%     fclose(fileID);
     
end

xlabel('\bf K (wave number)','FontSize',15)
ylabel('\bf \Omega_{imag}','FontSize',15)
xlim([0,(P+1)*pi])
xticks(0:pi:(P+1)*pi)
xticklabels({'\fontsize{14} 0','\fontsize{22}\pi'....
    ,'\fontsize{14}2\fontsize{22}\pi','\fontsize{14}3\fontsize{22}\pi'})

legend_cfl_max= strcat('CFL_{max}= ',num2str(CFL_max));
legend_cfl={'Semi-Discrete','0.1 CFL_{max}','0.2 CFL_{max}'....
    ,'0.3 CFL_{max}','0.4 CFL_{max}','0.5 CFL_{max}'....
    ,'0.6 CFL_{max}','0.7 CFL_{max}','0.8 CFL_{max}'.....
    ,'0.9 CFL_{max}',legend_cfl_max};

legend(legend_cfl,'FontSize',15)

dtitle = strcat('Fully-Discrete Dissipation analysis, mDGp'....
    ,num2str(P),' with  ',flux,', and RK',num2str(Prk).....
    ,'. \alpha=[',num2str(Alpha(1)),', '....
    ,num2str(Alpha(2)),', ',num2str(Alpha(3)), '], jump parameter');
title(dtitle,'FontSize',14);

set(gca,'Color',[0.8 0.8 0.8]);


%-------------------------------------------------------------------
%      fully disc Dispersion
%-------------------------------------------------------------------
h=figure;

for j=1:length(CFL)
    
     [~,DGfd]= mDG_FourStab(P,Prk, K, Beta, Alpha, CFL(j),true_tol);
     plot(K, DGfd.wp1(1,:),'color',ccolor_map(j,:),'linewidth',1.5),hold on
     
     % Saving the output:
%     fname = strcat(dir,'DGp',num2str(P),'_RK',num2str(Prk)....
%         ,'_Beta',num2str(Beta),'_CFL',num2str(CFL(j)),'_wp.dat');
%     fileID = fopen(fname,'w');
%     print_data= [K;DGfd.wp1(1,:);DGfd.wp2(1,:)];
%     fprintf(fileID,'%1.3f %1.10e %1.10e\r\n',print_data);
%     fclose(fileID);
end

plot(K,K,'--k')

xlabel('\bf K (wave number)','FontSize',15)
ylabel('\bf \Omega_{real}','FontSize',15)
xlim([0,(P+1)*pi])
xticks(0:pi:(P+1)*pi)
xticklabels({'\fontsize{14} 0','\fontsize{22}\pi'....
    ,'\fontsize{14}2\fontsize{22}\pi','\fontsize{14}3\fontsize{22}\pi'})

legend(legend_cfl,'FontSize',15)

ptitle = strcat('Fully-Discrete Dispersion analysis, mDGp'....
    ,num2str(P),' with  ',flux,', and RK',num2str(Prk).....
    ,'. \alpha=[',num2str(Alpha(1)),', '....
    ,num2str(Alpha(2)),', ',num2str(Alpha(3)), '], jump parameter');
title(ptitle,'FontSize',14);

set(gca,'Color',[0.8 0.8 0.8]);

clear figname

