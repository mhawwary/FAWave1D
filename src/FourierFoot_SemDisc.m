%======================================================================
% Function to compute the fourier foot print for the semi discrete Fourier 
% Stability analysis. 
% Fourier Foot Print consists of the Dispersion and Dissipation properties
% of a certain space discretization scheme. 
% This is only applied on the linear wave equation
%======================================================================

function [dissip,disper,L,V] = FourierFoot_SemDisc(A_sd)

    if(isscalar(A_sd))

        dissip = real(A_sd);
        disper = -imag(A_sd);

    else

        [V,E] = eig(A_sd);

        N = length(E(1,:));
        L = zeros(1,N);
        wr=L; wi=L;
        
        VV = V;

        for i=1:N

            L(1,i) = E(i,i);  % EigenValues

            wr(1,i) = real(L(1,i)); % Dissipation
            wi(1,i) = imag(L(1,i)); % Dispersion
            
            for j=1:N
                V(j,i) = V(j,i); 
            end
        end
        dissip = wr;
        disper = -wi;
    end

end