%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the Semi-Disc eqn for Compact Difference Schemes 
% Lele 1992, using modified wave-number analysis
% This is approximating the 1st derivative df/dx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Asd] = CD_semi_disc_eqn(order,Kwavenumber)
    K= Kwavenumber;

    if(order==4)   % 4th order compact scheme
        a = 3/2; b=0.0; alpha=1/4;
        Asd ....
            = -1i.*(a.*sin(K)+ (b/2).*sin(2.*K))./(1+(2.*alpha.*cos(K)));
        
    elseif(order==6)  % 6th order compact scheme
        a = 14/9; b=1/9; alpha=1/3;
        Asd ....
            = -1i.*(a.*sin(K)+ (b/2).*sin(2.*K))./(1+(2.*alpha.*cos(K)));
    end
end

