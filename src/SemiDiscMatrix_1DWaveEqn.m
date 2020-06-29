%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Mohammad Alhawwary, Uiversity of Kansas, 2017. 
% Function for DG Space Discretization scheme of the linear wave equation 
% This function gives you the semi-discrete amplification matrix 
% corressponding to DG scheme with polynomial orders up to p=8
% Order of Accuracy is p+1 
% Inputs: 
% Kwavenumber: wavenumber for exp(1i*K)
% Porder: polynomial order
% upwind_param: beta for hybrid numerical fluxes, i.e., beta \in [0,1],
% ref0: Cockburn&Shu2001,"Runge-Kutta Discontinuous Galerkin Methods for Convection-Dominated Problems",J.Sci.Comput.,doi:10.1023/A:1012873910884
% ref1: Alhawwary&Wang2018,"Fourier analysis and evaluation of DG, FD and Compact difference methods for conservation laws",j.comput.physics 2018, doi:10.1016/j.jcp.2018.07.018
% ref2: Alhawwary&Wang2018,  "Comparative Fourier Analysis of DG, FD and Compact Difference schemes", AIAA Aviation, 2018 Fluid Dynamics Conference, AIAA 2018-4267, doi:10.2514/6.2018-4267 
% Ouput: The semi-discrete matrix A, dUj/dt  = A Uj
function [Asd_DG] = SemiDiscMatrix_1DWaveEqn(Porder,Kwavenumber,upwind_param)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Beta = upwind_param; 
    K = Kwavenumber;
    nDOF= Porder+1;
    
    % These matrices was generated using Mathematica Notebooks
    A0 =[ 0, -1, 0, -1, 0, -1, 0, -1, 0;
        3, 0, -3, 0, -3, 0, -3, 0, -3;
        0, 5, 0, -5, 0, -5, 0, -5, 0;
        7, 0, 7, 0, -7, 0, -7, 0, -7;
        0, 9, 0, 9, 0, -9, 0, -9, 0;
        11, 0, 11, 0, 11, 0, -11, 0, -11;
        0, 13, 0, 13, 0, 13, 0, -13, 0;
        15, 0, 15, 0, 15, 0, 15, 0, -15;
        0, 17, 0, 17, 0, 17, 0, 17, 0];
    
    B0 =[-1, 0, -1, 0, -1, 0, -1, 0, -1;
        0, -3, 0, -3, 0, -3, 0, -3, 0;
        -5, 0, -5, 0, -5, 0, -5, 0, -5;
        0, -7, 0, -7, 0, -7, 0, -7, 0;
        -9, 0, -9, 0, -9, 0, -9, 0, -9;
        0, -11, 0, -11, 0, -11, 0, -11, 0;
        -13, 0, -13, 0, -13, 0, -13, 0, -13;
        0, -15, 0, -15, 0, -15, 0, -15, 0;
        -17, 0, -17, 0, -17, 0, -17, 0, -17 ];
    
    C0 = A0 + Beta.*B0 ;      % C0 * Uj
    
    Am1 =0.5.*[ 1, 1, 1, 1, 1, 1, 1, 1, 1;
        -3, -3, -3, -3, -3, -3, -3, -3, -3;
        5, 5, 5, 5, 5, 5, 5, 5, 5;
        -7, -7, -7, -7, -7, -7, -7, -7, -7;
        9, 9, 9, 9, 9, 9, 9, 9, 9;
        -11, -11, -11, -11, -11, -11, -11, -11, -11;
        13, 13, 13, 13, 13, 13, 13, 13, 13;
        -15, -15, -15, -15, -15, -15, -15, -15, -15;
        17, 17, 17, 17, 17, 17, 17, 17, 17];
    Bm1=0.5.* [ 1, 1, 1, 1, 1, 1, 1, 1, 1;
        -3, -3, -3, -3, -3, -3, -3, -3, -3;
        5, 5, 5, 5, 5, 5, 5, 5, 5;
        -7, -7, -7, -7, -7, -7, -7, -7, -7;
        9, 9, 9, 9, 9, 9, 9, 9, 9;
        -11, -11, -11, -11, -11, -11, -11, -11, -11;
        13, 13, 13, 13, 13, 13, 13, 13, 13;
        -15, -15, -15, -15, -15, -15, -15, -15, -15;
        17, 17, 17, 17, 17, 17, 17, 17, 17];
    
    Cm1 = Am1 + Beta.*Bm1;        % Cm1 * U_j-1

    Ap1= 0.5.*[ -1, 1, -1, 1, -1, 1, -1, 1, -1;
        -3, 3, -3, 3, -3, 3, -3, 3, -3;
        -5, 5, -5, 5, -5, 5, -5, 5, -5;
        -7, 7, -7, 7, -7, 7, -7, 7, -7;
        -9, 9, -9, 9, -9, 9, -9, 9, -9;
        -11, 11, -11, 11, -11, 11, -11, 11, -11;
        -13, 13, -13, 13, -13, 13, -13, 13, -13;
        -15, 15, -15, 15, -15, 15, -15, 15, -15;
        -17, 17, -17, 17, -17, 17, -17, 17, -17];
    Bp1= 0.5.* [ 1, -1, 1, -1, 1, -1, 1, -1, 1;
        3, -3, 3, -3, 3, -3, 3, -3, 3;
        5, -5, 5, -5, 5, -5, 5, -5, 5;
        7, -7, 7, -7, 7, -7, 7, -7, 7;
        9, -9, 9, -9, 9, -9, 9, -9, 9;
        11, -11, 11, -11, 11, -11, 11, -11, 11;
        13, -13, 13, -13, 13, -13, 13, -13, 13;
        15, -15, 15, -15, 15, -15, 15, -15, 15;
        17, -17, 17, -17, 17, -17, 17, -17, 17];
    
    Cp1 = Ap1 + Beta.* Bp1; % Cp1 * U_j+1

    DG_C0  = C0(1:nDOF,1:nDOF);
    DG_Cm1 = Cm1(1:nDOF,1:nDOF);
    DG_Cp1 = Cp1(1:nDOF,1:nDOF);
    
    Asd_DG =  (  DG_Cm1 .* exp(-1i.*K) + DG_C0 + DG_Cp1 .* exp(1i.*K) );

    % dU/dt = Cm1 * U_j-1 + C0 * U_j + Cp1 *U_j+1 ....(1)
    % and if periodic B.C. then: U_j+1 = exp(1i*K) * U_j, U_j-1 =
    % exp(-1i*K) * U_j, this assumes the intial condition is of the form
    % (U_j)_t=0 = exp(1i*K), and hence dU/dt= Asd_DG * Uj
    % If B.C. changes then you can use eqn(1) and substitute with the
    % appropeiate U_j+1, U_j-1 depending on the B.C., j: is the element
    % index, U_j=[u0,u1,...,up] DOFs nodal or modal. 
    
end