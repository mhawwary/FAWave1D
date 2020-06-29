clear
close all
clc
%--------------------------------------------------------------
% Script to plot DGp1  CFL for stability with different Beta's
%---------------------------------------------------------------
true_tol = 1.0;

P=1;           % Polynomial order
Prk = 2;       % RK order
bb = 0.0:0.05:1.0;
Nb = length(bb);
N=P+1;
Beta = ones(P+1,Nb);      % upwind_parameter

Beta(end,:) = bb;



CFL_test=0.005:0.005:1.00;
K_test = 0.005:0.005:(P+1)*pi;

[~,DGfdn]= mDG_FourStab(P,Prk, K_test, [1.00; 0.7], 0.445, 1.0);
[~,DGfd]= DG_FourStab(P,Prk, K_test, 0.80, 0.416,1.0); 
[~,DGfdu]=DG_FourStab(P,Prk, K_test, 1.00, 1./3.0,1.0);

figure
plot(K_test,DGfdn.wd1(1,:),'-k',K_test,DGfd.wd1(1,:)...
    ,'-r',K_test,DGfdu.wd1(1,:),'-.b')
grid on

title('\bf Dissipation ')
xlabel('K')
ylabel('\omega_{r}')
legend('DG, \beta_{1}=0.70, CFL=0.445','DG, all \beta=0.80, CFL=0.416','DG, all \beta=1.00, CFL=1/3')

figure
plot(K_test,DGfdn.wp1(1,:),'-k',K_test,DGfd.wp1(1,:)...
    ,'-r',K_test,DGfdu.wp1(1,:),'-.b'),hold on 
plot(K_test,K_test,'--k') 
grid on
title('\bf Dispersion ')
xlabel('K')
ylabel('\omega_{i}')
legend('DG, \beta_{1}=0.70, CFL=0.445','DG, all \beta=0.80, CFL=0.416','DG, all \beta=1.00, CFL=1/3')



CFL_max = zeros(1,Nb);
cfl_K = CFL_max;

for j=1:length(Beta)
    [CFL_max(j),cfl_K(j)] ....
        = identify_stable_CFL_DG_new....
        (P,Prk, K_test, Beta(:,j), CFL_test,1e-10);
end
 
% % Saving the output:
outdir='./matlab_results/DGp1/';

[max_cfl,ii] = max(CFL_max);

fprintf('Beta with Max CFL is, Beta:%1.2f, CFL:%1.3f\n',Beta(end,ii),max_cfl);

h=figure;

plot(Beta(end,:),CFL_max,'--or','linewidth',1.5)
xlabel('\bf\beta','FontSize',15)
ylabel('\bfCFL_{max}','FontSize',15)

dtitle = strcat('\bfStability limits for DGp'....
    ,num2str(P),' and RK',num2str(Prk),' vs upwind parameter \beta');
title(dtitle,'FontSize',14);

set(gca,'Color',[0.8 0.8 0.8]);




