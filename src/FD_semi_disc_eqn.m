%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the Semi-Disc eqn for Finite Difference Schemes 
% Using the explicit relation and Von-Neumann Analysis
% This is approximating the 1st derivative df/dx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Asd] = FD_semi_disc_eqn(order,Kwavenumber,upw_bias)

    K= Kwavenumber;
    
    if(order==1)  % 1st order upwind
        stencil_index=[0,-1];           %[j,j-1];
    elseif (order==2) % 2nd order central
        if(upw_bias==0)   % 2nd order central scheme
            stencil_index = [1,0,-1];       %[j+1,j,j-1];
        elseif(upw_bias==1)   % 2nd order fully-upwind
            stencil_index = [0,-1,-2]; 
        end
    elseif (order==3) % 3rd order
        if(upw_bias==0)   % fully upwind
            stencil_index = [0,-1,-2,-3];  %[j,j-1,j-2,j-3];
        elseif(upw_bias==1) % 1-point upwind biased
            stencil_index = [1,0,-1,-2];    %[j+1,j,j-1,j-2];
        end
    elseif (order==4) % 4th order central
        if(upw_bias==0)   % 4th order central scheme
            stencil_index = [2,1,0,-1,-2];  %[j+2,j+1,j,j-1,j-2];
        elseif(upw_bias==1)   % 2-point upwind-biased
            stencil_index = [1,0,-1,-2,-3];  %[j+2,j+1,j,j-1,j-2];
        elseif(upw_bias==2)   % fully-upwind
            stencil_index = [0,-1,-2,-3,-4];  %[j+2,j+1,j,j-1,j-2];
        end
    elseif (order==5) % 5th order
        if(upw_bias==0)   % fully upwind
            stencil_index = [0,-1,-2,-3,-4,-5];  %[j,j-1,...,j-5];
        elseif(upw_bias==1) % 3-point upwind biased
            stencil_index = [1,0,-1,-2,-3,-4];    %[j+1,j,j-1,...,j-4];
        elseif(upw_bias==2) % 1-point upwind biased
            stencil_index = [2,1,0,-1,-2,-3];    %[j+2,j+1,j,...,j-3];
        end
    elseif (order==6) % 6th order central
        if(upw_bias==0)   % 6th order central scheme
            stencil_index = [3,2,1,0,-1,-2,-3];  %[j+3,j+2,j+1,j,j-1,j-2,j-3];
        elseif(upw_bias==2) % 2-point upwind-biased
            stencil_index = [2,1,0,-1,-2,-3,-4];  %[j+3,j+2,j+1,j,j-1,j-2,j-3];
        elseif(upw_bias==3) % 4-point upwind-biased
            stencil_index = [1,0,-1,-2,-3,-4,-5];  %[j+3,j+2,j+1,j,j-1,j-2,j-3];
        elseif(upw_bias==1) % fully-upwind
            stencil_index = [0,-1,-2,-3,-4,-5,-6];  %[j+3,j+2,j+1,j,j-1,j-2,j-3];
        end
    elseif (order==7) % 7th order 
        if(upw_bias==0)   % fully upwind
            stencil_index = [0,-1,-2,-3,-4,-5,-6,-7];  %[j,j-1,...,j-5];
        elseif(upw_bias==1) % 3-point upwind biased
            stencil_index = [1,0,-1,-2,-3,-4,-5,-6];    %[j+1,j,j-1,...,j-4];
        elseif(upw_bias==2) % 2-point upwind biased
            stencil_index = [2,1,0,-1,-2,-3,-4,-5];    %[j+2,j+1,j,...,j-3];
        elseif(upw_bias==3) % 1-point upwind biased
            stencil_index = [3,2,1,0,-1,-2,-3,-4];    %[j+2,j+1,j,...,j-3];
        end
    end
    
    S = stencil_index;
    n=length(S);
    A = zeros(n,n);
    A(1,:)=ones(1,n);   
    B = zeros(n,1);
    B(2,1)=-1;

    for i=2:n
        for j=1:n
            A(i,j) = (- S(j))^(i-1) ; 
        end
    end

    FD_coeff =(A\B);
    C = zeros(1,n);
    for i=1:n
        C(1,i) = FD_coeff(i,1) .* exp(1i.*K.*S(i));
    end

    Asd = - sum(C);
      
end

