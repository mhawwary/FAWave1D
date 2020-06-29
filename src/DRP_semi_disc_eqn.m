%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the Semi-Disc eqn for DRP Schemes 
% Bogey and Bailly 2004, using modified wave-number analysis
% This is approximating the 1st derivative df/dx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Asd] = DRP_semi_disc_eqn(order,no_points,Kwavenumber)
    K= Kwavenumber;
    N = no_points;
    
    if(order==4)
        if(N==11)  % DRP(4,11)/FDo11p , 4th order of BogeyBailly
            a = [0.872756993962,-0.286511173973,0.090320001280....
                ,-0.020779405824,0.002484594688];
            summ = a(1).*sin(K)+a(2).*sin(2.*K)+a(3).*sin(3.*K)....
                              +a(4).*sin(4.*K)+a(5).*sin(5.*K); 
        elseif(N==13)  % DRP(4,13)/FDo13p , 4th order of BogeyBailly
            a = [0.907646591371,-0.337048393268,0.133442885327....
                ,-0.045246480208,0.011169294114,-0.001456501759];
            summ = a(1).*sin(K)+a(2).*sin(2.*K)+a(3).*sin(3.*K)....
                               +a(4).*sin(4.*K)+a(5).*sin(5.*K)....
                               +a(6).*sin(6.*K); 
        elseif(N==7)  % DRP(4,7)/FDo7p , 6th order of Tam
            a=[0.77088238051822552,-0.166705904414580469...
                ,0.02084314277031176 ]; % due to cunha
            a=[0.76868543,-0.16494834,0.02040375]; % due to Linders_Nordstrom2015
            summ = a(1).*sin(K)+a(2).*sin(2.*K)+a(3).*sin(3.*K);
        end
    elseif(order==2)  % Remez schemes introduced by Linders and Nordstrom
        if(N==7) % Remez(2,7)
            a=[0.78028389,-0.17585010,0.02380544];
            summ = a(1).*sin(K)+a(2).*sin(2.*K)+a(3).*sin(3.*K);
        elseif(N==9) % Remez(2,9)
            a=[0.82535056,-0.22669661,0.05059599,-0.00593632];
            summ = a(1).*sin(K)   +a(2).*sin(2.*K).....
                  +a(3).*sin(3.*K)+a(4).*sin(4.*K);
        end
    end
    
    Asd = -1i.*2.*summ;
    
    
end

