%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function is to define the basis functions for DG
function [alpha]=init_proj_coeff(P,phy,K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms Xi I

I = exp(1i.*K.*Xi/2).*phy ;
alpha = int(I,Xi,-1,1);
temp = int(phy.^2,Xi,-1,1);
alpha=alpha./temp;