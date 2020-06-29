function []= dump_dominant_mode(P,Prk,Beta,cfl_ratio,K,wd,wp)

if(Prk>0 || cfl_ratio>0)  % fully discrete
    outdir = strcat('../results/DGp',num2str(P),'_RK',num2str(Prk),'/');
    datadir='data/';
    fname = strcat(outdir,datadir,'DGp',num2str(P)....
        ,'_RK',num2str(Prk),'_Beta',num2str(Beta)...
        ,'_',num2str(cfl_ratio),'CFLmax_wdwp.dat');
else % semi-discrete
    outdir = strcat('../results/DGp',num2str(P),'/');
    datadir='data/';
    fname = strcat(outdir,datadir,'DGp',num2str(P)....
        ,'_Beta',num2str(Beta),'_sdisc_wdwp.dat');
end

fileID = fopen(fname,'w');
Nk = length(K);
print_data =zeros(Nk,3);

print_data(:,1)=K;
print_data(:,2)= wd;
print_data(:,3)= wp;

for j=1:Nk
    fprintf(fileID,'%1.3f %1.10e %1.10e\r\n',print_data(j,:));
end

fclose(fileID);
