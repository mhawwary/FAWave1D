%==========================================================================
% Function to compute the fourier Stability Properties of DG scheme and for
% different Polynomial order, numerical fluxes, and CFL no.
%==========================================================================

% upwind parameter = 1: upwind, 0: central
% Rk_order:  1: Euler, 2: RK2, 3: RK3, 4: RK4

function [DGsd,DGfd]....
    = mDG_FourStab(Porder,RK_order, Kwavenumber...
    , upwind_param, jump_param, CFL, true_tol)

K = Kwavenumber;
Beta=upwind_param;
Alpha=jump_param;
P=Porder;
Nk = length(K);
Nb = length(Beta(1,:));
Na = length(Alpha(1,:));

wd_name = cell(1,P+1);
wp_name = wd_name;

for i=1:P+1
    
    wd_name(1,i) = {strcat('wd', num2str(i))};
    wp_name(1,i) = {strcat('wp', num2str(i))};
    
end


wd_array = zeros(Nb,Nk);
wp_array = zeros(Nb,Nk);

for i=1:P+1
    DGsd.(char(wd_name(1,i)))=wd_array;
    DGsd.(char(wp_name(1,i)))=wp_array;
end

DGfd = DGsd;

 
%==========================================================================

for j=1:Nb   % Loop over Beta 
    for m=1:Na
        for k=1:Nk   % Loop over wave number

            %-----------------------------------------
            % Semi Discrete DG Stability Analysis
            %-----------------------------------------
            [Asd_DG] = mDG_semi_disc_matrix(P,K(k),Beta(:,j),Alpha(:,m)); 
            [temp_d, temp_p] = FourierFoot_SemDisc(Asd_DG);

            for i=1:P+1
                DGsd.(char(wd_name(1,i)))(j,k)=temp_d(i);
                DGsd.(char(wp_name(1,i)))(j,k)=temp_p(i);
            end

            %------------------------------
            % Fully Discrete DG Stability Analysis
            %------------------------------
            [Afd_DG]=RK(RK_order,Asd_DG,CFL);
            [temp_d, temp_p] = FourierFoot_FullDisc(Afd_DG,CFL);

            for i=1:P+1
                DGfd.(char(wd_name(1,i)))(j,k)=temp_d(i);
                DGfd.(char(wp_name(1,i)))(j,k)=temp_p(i);
            end

        end
        
        [DGsd] = DG_true_FourierFoot(DGsd,P,K,j,true_tol);
        [DGfd] = DG_true_FourierFoot(DGfd,P,K,j,true_tol);
        
    end
end

%==========================================================================

