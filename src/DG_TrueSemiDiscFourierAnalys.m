%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% True Semi-Discrete Fourier Analysis function
%
function [Kt,S,G, Kp,Sp,Gp, E, Eex]=....
    DG_TrueSemiDiscFourierAnalys(P,Beta,t,data_dump_flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nt = length(t);
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%                         True Semi Discrete Analysis                     %
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


if(data_dump_flag==0)
    ss=strcat('loading semi-disc DG','-Beta',num2str(Beta)...
        ,' scheme.........');  disp(ss);
    if(P==2)
%         Kp =[linspace(0.01,(P+1)*pi-0.1,50*(P+1)),(P+1)*pi]; % for true solution 
%         [DG,~]= DG_FourStab(P,3, Kp, Beta, 1, 2);
%         wd = DG.wd1(1,:)';
%         wp = DG.wp1(1,:)';
%         Kp=Kp';
        [Kp,wd,wp]= load_dominant_mode(P,0,Beta,0);
        figure
        plot(Kp./(P+1),wp./(P+1),'-r'),hold on
        
    elseif(P==5)
        [Kp,tempd,tempp]= load_dominant_mode_mat(P,0,Beta,0);
        
%         if(Beta==1)
%             tol = 3.0;
%         elseif(Beta==0)
%             tol=1.0;
%         end
%         Kp= [linspace(0.01,(P+1)*pi-0.1,50*(P+1)),(P+1)*pi]; % for true solution 
%         [DG,~]= DG_FourStab(P,3, Kp, Beta, 1, tol);
%         wd_name = cell(1,P+1);
%         wp_name = wd_name;
%         
%         figure
%         plot(Kp./(P+1),DG.wp1./(P+1),'-r'),hold on
%         plot(Kp./(P+1),DG.wp2./(P+1),'-k'),hold on
%         plot(Kp./(P+1),DG.wp3./(P+1),'-b'),hold on
%         plot(Kp./(P+1),DG.wp4./(P+1),'-g'),hold on
%         plot(Kp./(P+1),DG.wp5./(P+1),'-c'),hold on
%         plot(Kp./(P+1),DG.wp6./(P+1),'-m'),hold off
% 
%         for i=1:P+1
%             wd_name(1,i) = {strcat('wd', num2str(i))};
%             wp_name(1,i) = {strcat('wp', num2str(i))};   
%         end
%         for i=1:P+1
%            wd(:,i) = DG.(char(wd_name(1,i)))';
%            wp(:,i) = DG.(char(wp_name(1,i)))';
%         end
%         if(Beta==0)
%             tempd =  [wd(1,3);wd(2:92,4);wd(93:185,3);wd(186:284,1);wd(285:end,6)];
%             tempp =  [wp(1,3);wp(2:92,4);wp(93:185,3);wp(186:284,1);wp(285:end,6)];
%             
%         else
%             temp = wd(:,1);
%             wd(:,1) = wd(:,3);
%             wd(:,3) = temp;
%             temp = wp(:,1);
%             wp(:,1) = wp(:,3);
%             wp(:,3) = temp;
% 
%             temp = wd(110,1);
%             wd(110,1) = wd(110,6);
%             wd(110,6) = temp;
%             temp = wp(110,1);
%             wp(110,1) = wp(110,6);
%             wp(110,6) = temp;
%             temp = wd(111:117,1);
%             wd(111:117,1) = wd(111:117,5);
%             wd(111:117,5) = temp;
%             temp = wp(111:117,1);
%             wp(111:117,1) = wp(111:117,5);
%             wp(111:117,5) = temp;
% 
%             tempd=wd(:,1);
%             tempp=wp(:,1);
%         end
        figure
        plot(Kp./(P+1),tempp(:,1)./(P+1),'-r'),hold on
        wd = tempd;
        wp = tempp;
        %dump_SemiDisc_dominant_mode(P,Beta,Kp',wd,wp);
        %Kp=Kp';
    end
    Nk = length(Kp);
    Gp = zeros(Nt,Nk);
    Sp = Gp;
    
    for nt=1:length(t)
        ss=strcat('loading true semi-disc DG','-Beta',num2str(Beta)...
        ,' scheme.........');  disp(ss);
        [Kt,tempG,tempS]= load_true_G_S(P,3,Beta,0,t(nt));
        if(nt==1)
            G= zeros(Nt,length(Kt));
            S=G;
        end
        G(nt,:) = tempG;
        S(nt,:) = tempS;
        Gp(nt,:) = exp(wd.*t(nt));
        Sp(nt,:) = wp.*t(nt);
        disp(strcat('nt: ',num2str(nt)));
    end
    % need to be fixed
    E=0;
    Eex=0;
    
elseif(data_dump_flag==1)
    Kp=0;
    Sp=0;
    Gp=0;
    Kt= [linspace(0.01,(P+1)*pi-0.1,20*(P+1)),(P+1)*pi]; % for true solution 
    Nkt=length(Kt);
    Gv = zeros(Nkt,P+1); 
    phy=Gv; U = Gv; 
    g = zeros(P+1,P+1);
    G= zeros(Nt,Nkt);
    S = G;
    E=G; 
    Eex = G;
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
    
    for nt=1:length(t)
        [DGsd,~]= DG_FourStab(P,3, Kt, Beta, 1,0);
        for k=1:Nkt
            Ko=Kt(k);

            for i=1:P+1
                V(:,i) = DGsd.(char(V_name(1,i)))(:,k);
                L(i) =  DGsd.(char(L_name(1,i)))(1,k);
            end

            [alpha]=eval_init_proj_coeff(P,Ko);  % get the initial projection coefficients

            eta = V\alpha;  % modes weights
            I = 0;
            Ie =0;

            for j=1:N
                for i=1:N
                    g(j,i) = eta(i) * V(j,i) * exp(L(i)*t(nt));
                end
                U(k,j) = sum(g(j,:));
                Gv(k,j) = abs(U(k,j))/alpha(j,1);
                phy(k,j) = -angle(U(k,j));
                I = I+ U(k,j) * epsi(j);
                Ie = Ie + alpha(j,1) * epsi(j);
            end
            Ie = Ie * exp(-1i*Ko*t(nt));
            II = I * conj(Ie);
            I = abs(I)^2;
            Ie = abs(Ie)^2;
            E(nt,k) = sqrt(int(I,Xi,-1,1)/ 2);
            Eex(nt,k) = sqrt(int(Ie,Xi,-1,1)/ 2);
            II = int(II,Xi,-1,1);
            S(nt,k) = angle(II);
            G(nt,k) = E(nt,k)/Eex(nt,k);
        end
        disp(strcat('nt: ',num2str(nt)));
        dump_true_G_S(P,3,Beta,0,Kt,G(nt,:),S(nt,:),t(nt));
    end
end
