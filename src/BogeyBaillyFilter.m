function [A_bogey,Damp_func]=BogeyBaillyFilter(order,no_points,sigma_filter....
    ,Kwavenumber,G_fdisc)

K = Kwavenumber;
OA = order;
Np = no_points;
sigma=sigma_filter;

if(Np==11)
    d0 = 0.234810479761700;
    d= [-0.199250131285813,0.120198310245186.....
        ,-0.049303775636020,0.012396449873964,-0.001446093078167]; 
%     d0=0.215044884112;
%     d = [-0.187772883589,0.123755948787,-0.059227575576....
%         ,0.018721609157,-0.002999540835]; % old 2004 paper
    summ = d(1).*cos(K)+d(2).*cos(2.*K)+d(3).*cos(3.*K).....
                        +d(4).*cos(4.*K)+d(5).*cos(5.*K);
elseif(Np==13)
    d0 = 0.190899511506;
    d = [-0.171503832236,0.123632891797,-0.069975429105.....
        ,0.029662754736,-0.008520738659,0.001254597714];
    summ = d(1).*cos(K)+d(2).*cos(2.*K)+d(3).*cos(3.*K).....
                       +d(4).*cos(4.*K)+d(5).*cos(5.*K)+d(6).*cos(6.*K);
end

A_bogey = G_fdisc.*(1.0-sigma*d0-2.0*sigma.*summ);

Damp_func = d0 + 2.0*sigma.*summ;