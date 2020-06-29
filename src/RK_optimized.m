%======================================================================
% Function for Runge-Kutta time integration schemes for the 
% linear wave equation. It also includes the euler method. 
% This function gives you the amplification matrix corressponding to any
% space discretization scheme that it is coupled with RK. 
%======================================================================

function [G_fullyDiscrete] = RK_optimized(type,order,G_semiDiscrete,CFL)

    Asd =  G_semiDiscrete;

    if(strcmp(type,'CLASSICAL'))
        switch order
            case 1  % Euler forward method
                Afd = 1 +  CFL.* Asd;   
            case 2  % SSPRK(2,2)
                Afd = 1 +  CFL.* Asd ...
                    + ((CFL.^2)./2.0) .* (Asd * Asd);
            case 3  % SSPRK(3,3)
                Afd = 1 +  CFL.* Asd ...
                    + ((CFL.^2)./2.0) .* (Asd * Asd)...
                    + ((CFL.^3)./6.0) .* (Asd * Asd * Asd);
            case 4  % SSPRK(4,4)
                Afd = 1 +  CFL.* Asd ...
                    + ((CFL.^2)./2.0) .* (Asd * Asd)...
                    + ((CFL.^3)./6.0) .* (Asd * Asd * Asd)....
                    + ((CFL.^4)./24.0).* (Asd * Asd * Asd * Asd);
            otherwise
                disp('Wrong time integration order , or not implemented ')   
        end
        
    elseif(strcmp(type,'BogeyBailly6s'))  % RKo6s
        g = [1.00, 0.50, 0.165919771368,0.040919732041...
            ,0.007555704391,0.000891421261];
        Afd = 1 +g(1).*CFL.* Asd ...
                +g(2).*CFL.^2.* (Asd * Asd)...
                +g(3).*CFL.^3.* (Asd * Asd * Asd)....
                +g(4).*CFL.^4.* (Asd * Asd * Asd * Asd)....
                +g(5).*CFL.^5.* (Asd * Asd * Asd * Asd * Asd)....
                +g(6).*CFL.^6.* (Asd * Asd * Asd * Asd * Asd * Asd);
    end

    G_fullyDiscrete =Afd;

end