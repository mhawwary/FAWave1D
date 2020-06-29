%======================================================================
% Function to compute the fourier foot print for fully discrete Fourier 
% Stability analysis. 
% Fourier Foot Print consists of the Dispersion and Dissipation properties
% of a certain space and time schemes. 
% This is only applied on the linear wave equation
%======================================================================

function [dissip,disper] = FourierFoot_FullDisc(Lf,CFL)

    if(isscalar(Lf))

        D = (1i./CFL).*log(Lf);

        dissip = imag(D);
        disper = real(D);
        
        %disp('scalar in fourier foot fulldisc dg');

    else
        N = length(Lf);
        D =zeros(1,N);
        wr=D; wi = D;
        
        for i=1:N
            D(1,i) = (1i./CFL).* log(Lf(i));  % EigenValues
            wi(1,i) = imag(D(1,i)); % Dissipation
            wr(1,i) = real(D(1,i)); % Dispersion
        end

        dissip = wi;
        disper = wr;

    end


end