function [K,wd,wp]= load_dominant_mode(P,Prk,Beta,cfl_ratio)

if(Prk>0 || cfl_ratio>0)  % fully discrete
    indir = strcat('../results/DGp',num2str(P),'_RK',num2str(Prk),'/');
    datadir='data/';
    fname = strcat(indir,datadir,'DGp',num2str(P)....
        ,'_RK',num2str(Prk),'_Beta',num2str(Beta)...
        ,'_',num2str(cfl_ratio),'CFLmax_wdwp.dat');
else
    indir = strcat('../results/DGp',num2str(P),'/');
    datadir='data/';
    fname = strcat(indir,datadir,'DGp',num2str(P)....
        ,'_Beta',num2str(Beta),'_sdisc_wdwp.dat'); 
end

data=load(fname);
K=data(:,1);
wd=data(:,2);
wp=data(:,3);

