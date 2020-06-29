function [A_padee]=PadeeFilter(order,alpha_filter,Kwavenumber,G_fdisc)

K = Kwavenumber;
OA = order;
alpha=alpha_filter;

if(OA==6)
    a0 = (11/16) + (5/8)*alpha;
    a1 = (15/32) +(17/16)*alpha;
    a2 = (-3/16) +(3/8)*alpha;
    a3 = (1/32) - (1/16)*alpha;
    a4 = 0.0;
    a5 = 0.0;
    
elseif(OA==8)
    a0 = (93 + 70*alpha)/128;
    a1 = (7 + 18*alpha)/16;
    a2 = (-7 + 14*alpha)/32;
    a3 = (1/16) - (1/8)*alpha;
    a4 = (-1/128) +(1/64)*alpha;
    a5 = 0.0;
    
elseif(OA==10)
    a0 = (193 + 126*alpha)/256;
    a1 = (105 + 302*alpha)/256;
    a2 = 15*(-1 + 2*alpha)/64;
    a3 = 45*(1 - 2*alpha)/512;
    a4 = 5*(-1 + 2*alpha)/256;
    a5 = (1 - 2*alpha)/512;
    
end

A_padee = G_fdisc.*( a0 + a1.*cos(K) + a2.*cos(2.*K) + a3.*cos(3.*K) .....
    + a4.* cos(4.*K) + a5.* cos(5.*K) )....
        ./(1 + 2.* alpha.*cos(K) );