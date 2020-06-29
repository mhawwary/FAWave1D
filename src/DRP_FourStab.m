%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computez the fourier Stability Properties of DRP schemes
% for different Order of accuracy, and CFL no. 
% Both Semi- and Fully-Disc analysis are obtained
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% upwind parameter = 1: upwind, 0: central
% Rk_order:  1: Euler, 2: RK2, 3: RK3, 4: RK4  Classical schemes
% 5:RKo5s, 6:RKo6s  optimized schemes due to Bogey and Bailly 2004

function [SD,FD]= DRP_FourStab(OrderAccuracy, no_points.....
                               ,RK_type,RK_order,Kwavenumber,CFL)  

K = Kwavenumber;
Np = no_points;   % no. of stencil points
OA=OrderAccuracy;    % order of accuracy
Nk = length(K);

wd_array = zeros(1,Nk);
wp_array = zeros(1,Nk);

SD=struct('wd',wd_array,'wp',wp_array);  % semi-disc
FD=SD;   % fully-disc

%==========================================================================
for k=1:Nk   % Loop over wave number

    % Semi Discrete FD Stability Analysis
    %------------------------------
    [Asd] = DRP_semi_disc_eqn(OA, Np ,K(k));
    [SD.wd(1,k),SD.wp(1,k)] = FourierFoot_SemDisc(Asd); 
    
    if(CFL>0)
        % Fully Discrete FD Stability Analysis
        %------------------------------
        [Afd] = RK_optimized(RK_type, RK_order,Asd,CFL);
        [FD.wd(1,k),FD.wp(1,k)] = FourierFoot_FullDisc(Afd,CFL);
    end

end
%==========================================================================