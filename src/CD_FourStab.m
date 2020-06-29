%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computez the fourier Stability Properties of CD schemes
% for different Order of accuracy, and CFL no. 
% Both Semi- and Fully-Disc analysis are obtained
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% upwind parameter = 1: upwind, 0: central
% Rk_order:  1: Euler, 2: RK2, 3: RK3, 4: RK4

function [SD,FD]= CD_FourStab(OrderAccuracy....
    ,RK_order, Kwavenumber, CFL)  

K = Kwavenumber;
OA=OrderAccuracy;    % order of accuracy
Nk = length(K);

wd_array = zeros(1,Nk);
wp_array = zeros(1,Nk);

SD=struct('wd',wd_array,'wp',wp_array);  % semi-disc
FD=SD;   % fully-disc

%==========================================================================
for k=1:Nk   % Loop over wave number

    % Semi Discrete CD Stability Analysis
    %------------------------------
    [Asd] = CD_semi_disc_eqn(OA,K(k));
    [SD.wd(1,k),SD.wp(1,k)] = FourierFoot_SemDisc(Asd); 

    % Fully Discrete CD Stability Analysis
    %------------------------------
    [Afd] = RK(RK_order,Asd,CFL);
    [FD.wd(1,k),FD.wp(1,k)] = FourierFoot_FullDisc(Afd,CFL);

end
%==========================================================================
