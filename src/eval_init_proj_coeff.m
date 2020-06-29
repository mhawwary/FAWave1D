%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function is to evaluate the coefficients of the intial
%  solution projection of the fourier mode exp(1i * K) onto the Orthogonal
%  Legendre basis functions. If u= sum( U* phy)
function [alpha]=eval_init_proj_coeff(P,K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha=zeros(P+1,1);
z = K/2;

beta = (sin(z)./z) - cos(z);
alpha(1,1) = sqrt(2).*sin(z)./z; 
alpha(2,1) = (1i.*sqrt(6)./z).*beta;

for n=1:P-1
    a = sqrt(4*n+6)./z;
    b = mod(n,2);
    c = mod(n+1,2);
    I=0;
    for m=1:n
        d = mod(m+n+1,2);
        e = sqrt(m+0.5);
        I = I + e .* d .* alpha(m+1,1);
    end
    I = 1i.*I;
    alpha(n+2,1) = a .* (b .*sin(z) + 1i.*beta.*c + I);
end


temp = [2,2/3,2/5,2/7,2/9,2/11];  % phy norm2_squared
phy_norm2 = sqrt(temp);

for i=1:P+1
   alpha(i,1) = alpha(i,1)./phy_norm2(i); 
end