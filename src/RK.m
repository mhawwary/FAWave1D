%======================================================================
% Function for Runge-Kutta time integration schemes for the 
% linear wave equation. It also includes the euler method. 
% This function gives you the amplification matrix corressponding to any
% space discretization scheme that it is coupled with RK. 
%======================================================================

function [G_fullyDiscrete] = RK(order,G_semiDiscrete,CFL)

    Asd =  G_semiDiscrete;
    [n,m] = size(Asd);
    
    if(n==1 || m==1)  % eigenvalue modifications
        switch order
            case 1  % Euler forward method
                Afd = 1 +  CFL.* Asd; 
            case 2  % SSPRK(2,2)
                Afd = 1 +  CFL.* Asd ...
                    + ((CFL.^2)./2.0) .* (Asd .* Asd);
            case 3  % SSPRK(3,3)
                Afd = 1 +  CFL.* Asd...
                    + ((CFL.^2)./2.0) .* (Asd.* Asd)...
                    + ((CFL.^3)./6.0) .* (Asd.* Asd.* Asd);
            case 4  % SSPRK(4,4)
                Afd = 1 +  CFL.* Asd ...
                    + ((CFL.^2)./2.0) .* (Asd.* Asd)...
                    + ((CFL.^3)./6.0) .* (Asd.* Asd.* Asd)....
                    + ((CFL.^4)./24.0).* (Asd.* Asd.* Asd.* Asd);
            otherwise
                disp('Wrong time integration order , or not implemented ')   
        end
        
    else
        switch order
            case 1  % Euler forward method
                Afd = eye(size(Asd)) +  CFL.* Asd;   
            case 2  % SSPRK(2,2)
                Afd = eye(size(Asd)) +  CFL.* Asd ...
                    + ((CFL.^2)./2.0) .* (Asd * Asd);
            case 3  % SSPRK(3,3)
                Afd = eye(size(Asd)) +  CFL.* Asd ...
                    + ((CFL.^2)./2.0) .* (Asd * Asd)...
                    + ((CFL.^3)./6.0) .* (Asd * Asd * Asd);
            case 4  % SSPRK(4,4)
                Afd = eye(size(Asd)) +  CFL.* Asd ...
                    + ((CFL.^2)./2.0) .* (Asd * Asd)...
                    + ((CFL.^3)./6.0) .* (Asd * Asd * Asd)....
                    + ((CFL.^4)./24.0).* (Asd * Asd * Asd * Asd);
            otherwise
                disp('Wrong time integration order , or not implemented ')   
        end
    end
    
    G_fullyDiscrete =Afd;

end