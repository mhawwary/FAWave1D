%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% True Fully Discrete Fourier Analysis function
%
function [Kt,S,G, Kp,Sp,Gp, E, Eex]=....
    DG_TrueFullDiscFourierAnalys(P,Prk,Beta,CFL_max,ratio,Niter...
    ,data_dump_flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nt= length(Niter);
CFL = ratio .* CFL_max;

if(data_dump_flag==0)
    ss=strcat('loading DGp',num2str(P),'RK',num2str(Prk)...
        ,'-Beta',num2str(Beta),' scheme.........'); disp(ss);
    [Kp,wd,wp]=load_dominant_mode(P,Prk,Beta,ratio); 
    Nk = length(Kp);
    Gp = zeros(Nt,Nk);
    Sp = Gp;
%     G= zeros(Nt,10);
%     S=G;
    for nt=1:Nt
        ss=strcat('loading true semi-disc DG','-Beta',num2str(Beta)...
        ,' scheme.........');  disp(ss);
        [Kt,tempG,tempS]= load_true_G_S(P,Prk,Beta,ratio,Niter(nt));
        if(nt==1)
            G= zeros(Nt,length(Kt));
            S=G;
        end
        G(nt,:) = tempG;
        S(nt,:) = tempS;
        Gp(nt,:) = exp(wd.*Niter(nt).*CFL);
        Sp(nt,:) = wp.*Niter(nt)*CFL;
        disp(strcat('nt: ',num2str(nt)));
    end
    figure
    plot(Kp./(P+1),wp./(P+1),'-r')
    % need to be fixed
    E=0;
    Eex=0; 
    
elseif(data_dump_flag==1)
    
    Kt = [linspace(0.01,(P+1)*pi-0.1,20*(P+1)),(P+1)*pi]; % for true solution
    Nkt = length(Kt);
    Gv = zeros(Nkt,P+1); 
    phy=Gv; U = Gv; 
    g = zeros(P+1,P+1);
    G= zeros(Nt,Nkt);
    E=zeros(Nt,Nkt); 
    Eex = E;
    S = E;
    
    syms Xi epsi0 epsi1 epsi2 epsi3 epsi4 epsi5 epsi I Ie II
    epsi0 = Xi^(0);
    epsi1 = Xi;
    epsi2 = (1/2)*(3*Xi^2-1);
    epsi3 = (1/2)*(5*Xi^3-3*Xi);
    epsi4 = (1/8)*(35*Xi^4-30*Xi^2+3);
    epsi5 = (1/8)*(63*Xi^5-70*Xi^3+15*Xi);
    epsi=[epsi0, epsi1, epsi2, epsi3, epsi4, epsi5];

    L_name = cell(1,P+1);
    V_name = cell(1,P+1);

    for i=1:P+1
        L_name(1,i) = {strcat('L', num2str(i))};
        V_name(1,i) = {strcat('V', num2str(i))}; 
    end
    V = zeros(P+1,P+1);
    L = zeros(1,P+1);
    N=P+1;

    for nt=1:Nt   % niter loop
        [~,DGfd]= DG_FourStab(P,Prk, Kt, Beta, CFL,0);
        for k=1:Nkt  % Kt loop
            %disp(strcat('n: ',num2str(n)));
            Ko=Kt(k);
            [alpha]=eval_init_proj_coeff(P,Ko);

            for i=1:P+1
                V(:,i) = DGfd.(char(V_name(1,i)))(:,k);
                L(i) =  DGfd.(char(L_name(1,i)))(1,k);
            end

            eta = V\alpha;  % modes weights
            I = 0;
            Ie =0;

            for j=1:N
                for i=1:N
                    g(j,i) = eta(i) * V(j,i) * L(i)^Niter(nt);
%                     Km(k,i) = (1i./CFL).* log(L(i));  % Km 
%                     wi(k,i) = imag(Km(k,i)); % Dissipation
%                     wr(k,i) = real(Km(k,i)); % Dispersion
                end
                U(k,j) = sum(g(j,:));
                Gv(k,j) = abs(U(k,j))/abs(alpha(j,1));
                phy(k,j) = -angle(U(k,j));
                I = I+ U(k,j) * epsi(j);
                Ie = Ie + alpha(j,1) * epsi(j);
            end
            Ie = Ie .* exp(-1i*Ko*CFL*Niter(nt));
            II = I .* conj(Ie);
            I = abs(I).^2;
            Ie = abs(Ie).^2;
            I = int(I,Xi,-1,1);
            E(nt,k) = sqrt(I .* 0.5);
            Ie = int(Ie,Xi,-1,1);
            Eex(nt,k) = sqrt(Ie .* 0.5);
            II = int(II,Xi,-1,1);
            S(nt,k) = angle(II);
            G(nt,k) = E(nt,k)/Eex(nt,k);
        end
        disp(strcat('nt: ',num2str(nt)));
        dump_true_G_S(P,Prk,Beta,ratio,Kt,G(nt,:),S(nt,:),Niter(nt));
    end
    
    Kp=0.0;
    Gp=0.0;
    Sp=0.0;

end

