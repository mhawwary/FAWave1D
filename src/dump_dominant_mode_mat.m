function []= dump_dominant_mode_mat(P,Prk,Beta,cfl_ratio,Kp,wdp,wpp)
K=Kp;
wd=wdp;
wp=wpp;

if(Prk>0 || cfl_ratio>0)  % fully discrete
    outdir = strcat('../results/DGp',num2str(P),'_RK',num2str(Prk),'/');
    datadir='data/';
    fname = strcat(outdir,datadir,'DGp',num2str(P)....
        ,'_RK',num2str(Prk),'_Beta',num2str(Beta)...
        ,'_',num2str(cfl_ratio),'CFLmax_wdwp.mat');
else % semi-discrete
    outdir = strcat('../results/DGp',num2str(P),'/');
    datadir='data/';
    fname = strcat(outdir,datadir,'DGp',num2str(P)....
        ,'_Beta',num2str(Beta),'_sdisc_wdwp.mat');
end
save(fname,'K','wd','wp','-mat');
