%==========================================================================
% Function to compute the fourier Stability Properties of DG scheme and for
% different Polynomial order, numerical fluxes, and CFL no.
%==========================================================================

% upwind parameter = 1: upwind, 0: central
% Rk_order:  1: Euler, 2: RK2, 3: RK3, 4: RK4

% Until now this function will not allow for testing eigenvectors of
% different Beta's, only one Beta is allowed per call for this purpose
% , however, if one needs to just steady L, wd, wp then it is possible

function [DGsd,DGfd]....
    = DG_FourStab(Porder,RK_order, Kwavenumber...
    , upwind_param, CFL, true_tol)

K = Kwavenumber;
Beta=upwind_param;
P=Porder;
Nk = length(K);
Nb = length(Beta);
Nm =P+1;  % no. of modes

wd_name = cell(1,Nm);
wp_name = wd_name;
L_name = wd_name;
V_name = cell(1,Nm);
we_name = wd_name;
alpha_name=wd_name;

for i=1:Nm
    wd_name(1,i) = {strcat('wd', num2str(i))};
    wp_name(1,i) = {strcat('wp', num2str(i))};
    L_name(1,i) = {strcat('L', num2str(i))};
    V_name(1,i) = {strcat('V', num2str(i))}; 
    we_name(1,i) = {strcat('we', num2str(i))};
    alpha_name(1,i) = {strcat('alpha', num2str(i))};
end

wd_array = zeros(Nb,Nk);
V_array = zeros(Nm,Nk);

for i=1:Nm
    DGsd.(char(wd_name(1,i)))=wd_array;
    DGsd.(char(wp_name(1,i)))=wd_array;
    DGsd.(char(L_name(1,i)))=wd_array;
    DGsd.(char(V_name(1,i)))=V_array;
    DGsd.(char(we_name(1,i)))=wd_array;
    DGsd.(char(alpha_name(1,i)))=wd_array;
end
DGfd = DGsd;
%==========================================================================

for j=1:Nb   % Loop over Beta 
    for k=1:Nk   % Loop over wavenumber K
        %------------------------------
        % Semi Discrete DG Stability Analysis
        %------------------------------
        [Asd] = SemiDiscMatrix_1DWaveEqn(P,K(k),Beta(j)); 
        [temp_d,temp_p,L,V] = FourierFoot_SemDisc(Asd);
        
        [alpha]=eval_init_proj_coeff(P,K(k));  % init_proj_coeff
        theta=V\alpha;
        for i=1:Nm
            DGsd.(char(wd_name(1,i)))(j,k)=temp_d(i);
            DGsd.(char(wp_name(1,i)))(j,k)=temp_p(i);
            DGsd.(char(L_name(1,i)))(j,k)=L(i);
            for m = 1:Nm
                DGsd.(char(V_name(1,i)))(m,k)=V(m,i);
            end
            DGsd.(char(we_name(1,i)))(j,k) = theta(i);  % modes weights
            DGsd.(char(alpha_name(1,i)))(j,k) = alpha(i); % init_proj_coeff
        end

        %------------------------------
        % Fully Discrete DG Stability Analysis
        %------------------------------
        if(RK_order>0)
            [Lf] = RK(RK_order,L,CFL);
            [temp_d, temp_p] = FourierFoot_FullDisc(Lf,CFL);

            for i=1:Nm
                DGfd.(char(wd_name(1,i)))(j,k)=temp_d(i);
                DGfd.(char(wp_name(1,i)))(j,k)=temp_p(i);

                DGfd.(char(L_name(1,i)))(j,k)=Lf(i);
                for m = 1:P+1
                    DGfd.(char(V_name(1,i)))(m,k)=V(m,i);
                end
                DGfd.(char(we_name(1,i)))(j,k)....
                    =DGsd.(char(we_name(1,i)))(j,k);
                DGfd.(char(alpha_name(1,i)))(j,k)....
                    =DGsd.(char(alpha_name(1,i)))(j,k);
            end
            
        end

    end
    
    if(true_tol~=0)
        if(RK_order>0)
            [DGfd] = DG_PhysFourierFoot(DGfd,P,K,j,true_tol);
        else
            [DGsd] = DG_PhysFourierFoot(DGsd,P,K,j,true_tol);
        end
    end
    
end

%==========================================================================

