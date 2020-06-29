function []= dump_true_G_S(P,Prk,Beta,cfl_ratio,K,G,S,time)


datadir='data/';

if(cfl_ratio==0)
    outdir = strcat('../results/DGp',num2str(P),'/');
    fname = strcat(outdir,datadir,'p',num2str(P)....
        ,'_Beta',num2str(Beta)...
        ,'_sdisc_G_S_t',num2str(time),'.dat');
elseif(cfl_ratio>0)
    outdir = strcat('../results/DGp',num2str(P),'_RK',num2str(Prk),'/');
    fname = strcat(outdir,datadir,'p',num2str(P)....
        ,'RK',num2str(Prk),'_Beta',num2str(Beta)...
        ,'_',num2str(cfl_ratio),'CFLmax_G_S_iter',num2str(time),'.dat');
end

fileID = fopen(fname,'w');
Nk = length(K);
print_data =zeros(Nk,3);

print_data(:,1)=K;
print_data(:,2)= G;
print_data(:,3)= S;

for j=1:Nk
    fprintf(fileID,'%1.3f %1.10e %1.10e\r\n',print_data(j,:));
end

fclose(fileID);
