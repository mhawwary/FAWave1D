function [K,G,S]= load_true_G_S(P,Prk,Beta,cfl_ratio,time)

datadir='data/';

if(cfl_ratio==0)
    indir = strcat('../results/DGp',num2str(P),'/');
    fname = strcat(indir,datadir,'p',num2str(P)....
        ,'_Beta',num2str(Beta)...
        ,'_sdisc_G_S_t',num2str(time),'.dat');
elseif(cfl_ratio>0)
    indir = strcat('../results/DGp',num2str(P),'_RK',num2str(Prk),'/');
    fname = strcat(indir,datadir,'p',num2str(P)....
        ,'RK',num2str(Prk),'_Beta',num2str(Beta)...
        ,'_',num2str(cfl_ratio),'CFLmax_G_S_iter',num2str(time),'.dat');
end

data=load(fname);
K=data(:,1);
G=data(:,2);
S=data(:,3);
