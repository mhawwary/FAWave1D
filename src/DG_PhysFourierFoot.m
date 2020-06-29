%==========================================================================
% This Function helps to distinguish different modes from one other
% and be able to plot them in a nice way
%==========================================================================

function [DG] = DG_PhysFourierFoot(DG,Porder,K,Beta_index,tol)

% tol is in degree

    j= Beta_index;

    % Check jumps :

    P=Porder;

    wd_name = cell(1,P+1);
    wp_name = wd_name;
    we_name = wd_name;
    V_name = wd_name;

    for i=1:P+1
        wd_name(1,i) = {strcat('wd', num2str(i))};
        wp_name(1,i) = {strcat('wp', num2str(i))};
        we_name(1,i) = {strcat('we', num2str(i))};
        V_name(1,i) = {strcat('V', num2str(i))};
    end

    % calculating the initial slope for each curve:
    %-----------------------------------------------
    k=2;
    slope_old=zeros(1,P+1);
    
    for m=1:P+1
        dDG_old = DG.(char(wp_name(1,m)))(j,k+1)-DG.(char(wp_name(1,m)))(j,k-1); 
        slope_old(m) = atan2( dDG_old, 2*(K(k+1)-K(k-1)) ) * 180/pi;
    end
    
    for k=2:length(K)-1

        %check_jump=logical(false);   % check for jump
        %check_flip=logical(false);   % check if slope change sign apprubtly
        %check_change=logical(false); % check if slope change in a high amount 
        
        p=1;
        n=1;
        
        dK = K(k+1)-K(k);
        d2K = 2 * (K(k+1)-K(k-1));
        
        while (p<P+1)

            dDG = DG.(char(wp_name(1,p)))(j,k+1)....
                -DG.(char(wp_name(1,p)))(j,k); 
            
            slope = atan2( dDG,dK ) * 180/pi;
            
            check_error = abs( slope-slope_old(p) ) > 2.5; % 2.5
            check_sign = ( slope*slope_old(p)) < 0;
            
            if(check_sign) 
                
                dDG = DG.(char(wp_name(1,p)))(j,k+1)....
                -DG.(char(wp_name(1,p)))(j,k-1);
                
                slope_2 = atan2( dDG,d2K ) * 180/pi;
                
                check_extrema = abs(slope_2) < 0.5;
                
                if(check_extrema)  % it is an extrema
                    check_sign = false;
                    check_error = false;
                else
                    check_sign = true;
                    check_error = true;
                end
            end
            
            if ( check_error || check_sign)  %  if there is a jump
                
                for ip=p+1:P+1

                    dDG = DG.(char(wp_name(1,ip)))(j,k+1).....
                        -DG.(char(wp_name(1,p)))(j,k); 
                    
                    slope=atan2( dDG, dK ) * 180/pi;
                    
                    check_error = abs( slope-slope_old(p) ) > 2.5;
                    check_sign = ( slope*slope_old(p)) < 0;

                    if(check_sign) 
                        
                        dDG = DG.(char(wp_name(1,ip)))(j,k+1)....
                        -DG.(char(wp_name(1,p)))(j,k-1);

                        slope_2 = atan2( dDG,d2K ) * 180/pi;

                        check_extrema = abs(slope_2) < 0.5;
                        
                         if(check_extrema)  % it is an extrema
                            check_sign = false;
                            check_error = false;
                        else
                            check_sign = true;
                            check_error = true;
                        end
                    end
                    
                    if ( check_error || check_sign)  % true if there is a jump

                    else
                        temp_old = DG.(char(wp_name(1,p)))(j,k+1:end);
                        temp_new = DG.(char(wp_name(1,ip)))(j,k+1:end);

                        DG.(char(wp_name(1,p)))(j,k+1:end) = temp_new;
                        DG.(char(wp_name(1,ip)))(j,k+1:end) = temp_old;

                        temp_old = DG.(char(wd_name(1,p)))(j,k+1:end);
                        temp_new = DG.(char(wd_name(1,ip)))(j,k+1:end);

                        DG.(char(wd_name(1,p)))(j,k+1:end) = temp_new;
                        DG.(char(wd_name(1,ip)))(j,k+1:end) = temp_old;
                        
                        temp_old = DG.(char(we_name(1,p)))(j,k+1:end);
                        temp_new = DG.(char(we_name(1,ip)))(j,k+1:end);

                        DG.(char(we_name(1,p)))(j,k+1:end) = temp_new;
                        DG.(char(we_name(1,ip)))(j,k+1:end) = temp_old;
                        
                        for ss=1:P+1
                            temp_old = DG.(char(V_name(1,p)))(ss,k+1:end);
                            temp_new = DG.(char(V_name(1,ip)))(ss,k+1:end);

                            DG.(char(V_name(1,p)))(ss,k+1:end) = temp_new;
                            DG.(char(V_name(1,ip)))(ss,k+1:end) = temp_old;
                        end

                        dDG = DG.(char(wp_name(1,ip)))(j,k+1).....
                            -DG.(char(wp_name(1,ip)))(j,k); 

                        slope_old(ip)= atan2( dDG, dK ) * 180/pi;

                        dDG = DG.(char(wp_name(1,p)))(j,k+1).....
                            -DG.(char(wp_name(1,p)))(j,k); 

                        slope_old(p)= atan2( dDG, dK ) * 180/pi;
                        
                        break;
                    end
                  
                end
                
            else   
                slope_old(p)= slope;
            end
            
            dDG = DG.(char(wp_name(1,p)))(j,k+1).....
                            -DG.(char(wp_name(1,p)))(j,k); 

            slope_old(p)= atan2( dDG, dK ) * 180/pi;
    
            p=p+1;
            n=n+1; 
        end
        
        dDG = DG.(char(wp_name(1,P+1)))(j,k+1).....
                            -DG.(char(wp_name(1,P+1)))(j,k); 

        slope_old(P+1)= atan2( dDG, dK ) * 180/pi;
        
        repeat=0;
        jump_mode = zeros(1,P+1);
        
        for ii=1:P+1
            
            if(abs(abs(slope_old(ii)) - 90) <=tol)
                
%                 disp(strcat('K: ',num2str(K(k)),', and p= ',num2str(ii)....
%                     ,', slope:',num2str(slope_old(ii)))) ;
                
                jump_mode(ii) = 1;
                
                repeat=repeat+1;
                
            end
            
        end
        
        ff= find(jump_mode==1);
        if(length(ff) > 2 )
            disp('problem')
        end
        
        if(repeat>1)
            
            p= ff(1);
            ip= ff(2);
            
            temp_old = DG.(char(wp_name(1,p)))(j,k+1:end);
            temp_new = DG.(char(wp_name(1,ip)))(j,k+1:end);

            DG.(char(wp_name(1,p)))(j,k+1:end) = temp_new;
            DG.(char(wp_name(1,ip)))(j,k+1:end) = temp_old;

            temp_old = DG.(char(wd_name(1,p)))(j,k+1:end);
            temp_new = DG.(char(wd_name(1,ip)))(j,k+1:end);

            DG.(char(wd_name(1,p)))(j,k+1:end) = temp_new;
            DG.(char(wd_name(1,ip)))(j,k+1:end) = temp_old;
            
            temp_old = DG.(char(we_name(1,p)))(j,k+1:end);
            temp_new = DG.(char(we_name(1,ip)))(j,k+1:end);

            DG.(char(we_name(1,p)))(j,k+1:end) = temp_new;
            DG.(char(we_name(1,ip)))(j,k+1:end) = temp_old;

            dDG = DG.(char(wp_name(1,ip)))(j,k+1).....
                -DG.(char(wp_name(1,ip)))(j,k); 

            slope_old(ip)= atan2( dDG, dK ) * 180/pi;

            dDG = DG.(char(wp_name(1,p)))(j,k+1).....
                -DG.(char(wp_name(1,p)))(j,k); 

            slope_old(p)= atan2( dDG, dK ) * 180/pi;
                 
        end

    end

    
    
end
