function [CFL_max,K_unstable_cfl]= identify_stable_CFL_DG(P,Prk, K, Beta, CFL,tol)

check1=0;
CFL_max=0;
K_unstable_cfl=0;

wd_name = cell(1,P+1);

for i=1:P+1
    wd_name(1,i) = {strcat('wd', num2str(i))};
end

for j=1:length(CFL)
    
    [~,DGfd]= DG_FourStab(P,Prk, K, Beta, CFL(j),0.0);
    
    i=1;
    
    if(check1==0)
        
        while( (check1==0 ) && i<length(K) )
            for ii=1:P+1
                
                wd = DGfd.(char(wd_name(1,ii)))(1,i);
                
                if(wd > tol ) 
                    
                    if(j==1) 
                        CFL_max=CFL(j);
                    else
                        CFL_max = CFL(j-1);
                    end
                    pp = strcat('DGp',num2str(P),' and RK',num2str(Prk)....
                        ,' Beta= ',num2str(Beta));
                    disp(pp);
                    pCFL = strcat('Unstable CFL: ',num2str(CFL(j))...
                        ,', at K: ',num2str(K(i))...
                        ,', |G|= ',num2str(exp(CFL(j)*wd))....
                        ,', Omega: ',num2str(wd));
                    disp(pCFL);
                    maxCFL = strcat('CFL_max: ',num2str(CFL_max));
                    disp(maxCFL);
                    K_unstable_cfl=K(i);
                    check1=1; 
                    break;
                end
            end
            i=i+1;
        end
            
%             i=1;
%             while( (check1==0 ) && i<length(K) )
%                 if(DGfd.wd2(1,i) > tol ) 
%                     if(j==1) 
%                         CFL_max=CFL(j);
%                     else
%                         CFL_max = CFL(j-1);
%                     end
%                     
%                     pp = strcat('DGp',num2str(P),' and RK',num2str(Prk),' Beta= ',num2str(Beta));
%                     disp(pp);
%                     pCFL = strcat('Unstable CFL: ',num2str(CFL(j))...
%                         ,', at K: ',num2str(K(i))...
%                         ,', |G|= ',num2str(exp(CFL(j)*DGfd.wd2(1,i))),', Omega1: ',num2str(DGfd.wd2(1,i)));
%                     disp(pCFL);
%                     maxCFL = strcat('CFL_max: ',num2str(CFL_max));
%                     disp(maxCFL);
%                     K_unstable_cfl=K(i);
%                     check1=1; 
%                 end
%                 i=i+1;
%              end
% 
%         elseif(P==2)
%             
%             while( (check1==0 ) && i<length(K) )
%                 if(DGfd.wd1(1,i) > tol ) 
%                     if(j==1) 
%                         CFL_max=CFL(j);
%                     else
%                         CFL_max = CFL(j-1);
%                     end
%                     pp = strcat('DGp',num2str(P),' and RK',num2str(Prk),' Beta= ',num2str(Beta));
%                     disp(pp);
%                     pCFL = strcat('Unstable CFL: ',num2str(CFL(j))...
%                         ,', at K: ',num2str(K(i))...
%                         ,', |G|= ',num2str(exp(CFL(j)*DGfd.wd1(1,i))),', Omega1: ',num2str(DGfd.wd1(1,i)));
%                     disp(pCFL);
%                     maxCFL = strcat('CFL_max: ',num2str(CFL_max));
%                     disp(maxCFL);
%                     K_unstable_cfl=K(i);
%                     check1=1; 
%                 end
%                 i=i+1;
%             end
%             
%             i=1;
%             while( (check1==0 ) && i<length(K) )
%                 if(DGfd.wd2(1,i) > tol ) 
%                     if(j==1) 
%                         CFL_max=CFL(j);
%                     else
%                         CFL_max = CFL(j-1);
%                     end
%                     pp = strcat('DGp',num2str(P),' and RK',num2str(Prk),' Beta= ',num2str(Beta));
%                     disp(pp);
%                     pCFL = strcat('Unstable CFL: ',num2str(CFL(j))...
%                         ,', at K: ',num2str(K(i))...
%                         ,', |G|= ',num2str(exp(CFL(j)*DGfd.wd2(1,i))),', Omega2: ',num2str(DGfd.wd2(1,i)));
%                     disp(pCFL);
%                     maxCFL = strcat('CFL_max: ',num2str(CFL_max));
%                     disp(maxCFL);
%                     K_unstable_cfl=K(i);
%                     check1=1; 
%                 end
%                 i=i+1;
%             end
%             
%             i=1;
%             while( (check1==0 ) && i<length(K) )
%                 if(DGfd.wd3(1,i) > tol ) 
%                     if(j==1) 
%                         CFL_max=CFL(j);
%                     else
%                         CFL_max = CFL(j-1);
%                     end
%                     pp = strcat('DGp',num2str(P),' and RK',num2str(Prk),' Beta= ',num2str(Beta));
%                     disp(pp);
%                     pCFL = strcat('Unstable CFL: ',num2str(CFL(j))...
%                         ,', at K: ',num2str(K(i))...
%                         ,', |G|= ',num2str(exp(CFL(j)*DGfd.wd3(1,i))),', Omega3: ',num2str(DGfd.wd3(1,i)));
%                     disp(pCFL);
%                     maxCFL = strcat('CFL_max: ',num2str(CFL_max));
%                     disp(maxCFL);
%                     K_unstable_cfl=K(i);
%                     check1=1; 
%                 end
%                 i=i+1;
%             end 
% 
%         elseif(P==3)
%             
%             while( (check1==0 ) && i<length(K) )
%                 if(DGfd.wd1(1,i) > tol ) 
%                      if(j==1) 
%                         CFL_max=CFL(j);
%                     else
%                         CFL_max = CFL(j-1);
%                     end
%                     pp = strcat('DGp',num2str(P),' and RK',num2str(Prk),' Beta= ',num2str(Beta));
%                     disp(pp);
%                     pCFL = strcat('Unstable CFL: ',num2str(CFL(j))...
%                         ,', at K: ',num2str(K(i))...
%                         ,', |G|= ',num2str(exp(CFL(j)*DGfd.wd1(1,i))),', Omega1: ',num2str(DGfd.wd1(1,i)));
%                     disp(pCFL);
%                     maxCFL = strcat('CFL_max: ',num2str(CFL_max));
%                     disp(maxCFL);
%                     K_unstable_cfl=K(i);
%                     check1=1; 
%                 end
%                 i=i+1;
%             end
%              
%             while( (check1==0 ) && i<length(K) )
%                 if(DGfd.wd2(1,i) > tol ) 
%                      if(j==1) 
%                         CFL_max=CFL(j);
%                     else
%                         CFL_max = CFL(j-1);
%                     end
%                     pp = strcat('DGp',num2str(P),' and RK',num2str(Prk),' Beta= ',num2str(Beta));
%                     disp(pp);
%                     pCFL = strcat('Unstable CFL: ',num2str(CFL(j))...
%                         ,', at K: ',num2str(K(i))...
%                         ,', |G|= ',num2str(exp(CFL(j)*DGfd.wd2(1,i))),', Omega1: ',num2str(DGfd.wd2(1,i)));
%                     disp(pCFL);
%                     maxCFL = strcat('CFL_max: ',num2str(CFL_max));
%                     disp(maxCFL);
%                     K_unstable_cfl=K(i);
%                     check1=1; 
%                 end
%                 i=i+1;
%             end
%              
%             while( (check1==0 ) && i<length(K) )
%                 if(DGfd.wd3(1,i) > tol ) 
%                      if(j==1) 
%                         CFL_max=CFL(j);
%                     else
%                         CFL_max = CFL(j-1);
%                     end
%                     pp = strcat('DGp',num2str(P),' and RK',num2str(Prk),' Beta= ',num2str(Beta));
%                     disp(pp);
%                     pCFL = strcat('Unstable CFL: ',num2str(CFL(j))...
%                         ,', at K: ',num2str(K(i))...
%                         ,', |G|= ',num2str(exp(CFL(j)*DGfd.wd3(1,i))),', Omega1: ',num2str(DGfd.wd3(1,i)));
%                     disp(pCFL);
%                     maxCFL = strcat('CFL_max: ',num2str(CFL_max));
%                     disp(maxCFL);
%                     K_unstable_cfl=K(i);
%                     check1=1; 
%                 end
%                 i=i+1;
%             end 
%             
%             while( (check1==0 ) && i<length(K) )
%                 if(DGfd.wd4(1,i) > tol ) 
%                      if(j==1) 
%                         CFL_max=CFL(j);
%                     else
%                         CFL_max = CFL(j-1);
%                     end
%                     pp = strcat('DGp',num2str(P),' and RK',num2str(Prk),' Beta= ',num2str(Beta));
%                     disp(pp);
%                     pCFL = strcat('Unstable CFL: ',num2str(CFL(j))...
%                         ,', at K: ',num2str(K(i))...
%                         ,', |G|= ',num2str(exp(CFL(j)*DGfd.wd4(1,i))),', Omega1: ',num2str(DGfd.wd4(1,i)));
%                     disp(pCFL);
%                     maxCFL = strcat('CFL_max: ',num2str(CFL_max));
%                     disp(maxCFL);
%                     K_unstable_cfl=K(i);
%                     check1=1; 
%                 end
%                 i=i+1;
%             end 
%             
%         elseif(P==4)
%             
%             while( (check1==0 ) && i<length(K) )
%                 if(DGfd.wd1(1,i) > tol ) 
%                     CFL_max = CFL(j-1);
%                     pp = strcat('DGp',num2str(P),' and RK',num2str(Prk),' Beta= ',num2str(Beta));
%                     disp(pp);
%                     pCFL = strcat('Unstable CFL: ',num2str(CFL(j))...
%                         ,', at K: ',num2str(K(i))...
%                         ,', |G|= ',num2str(exp(CFL(j)*DGfd.wd1(1,i))),', Omega1: ',num2str(DGfd.wd1(1,i)));
%                     disp(pCFL);
%                     maxCFL = strcat('CFL_max: ',num2str(CFL_max));
%                     disp(maxCFL);
%                     K_unstable_cfl=K(i);
%                     check1=1; 
%                 end
%                 i=i+1;
%             end
%              
%             while( (check1==0 ) && i<length(K) )
%                 if(DGfd.wd2(1,i) > tol ) 
%                     CFL_max = CFL(j-1);
%                     pp = strcat('DGp',num2str(P),' and RK',num2str(Prk),' Beta= ',num2str(Beta));
%                     disp(pp);
%                     pCFL = strcat('Unstable CFL: ',num2str(CFL(j))...
%                         ,', at K: ',num2str(K(i))...
%                         ,', |G|= ',num2str(exp(CFL(j)*DGfd.wd2(1,i))),', Omega1: ',num2str(DGfd.wd2(1,i)));
%                     disp(pCFL);
%                     maxCFL = strcat('CFL_max: ',num2str(CFL_max));
%                     disp(maxCFL);
%                     K_unstable_cfl=K(i);
%                     check1=1; 
%                 end
%                 i=i+1;
%             end
%              
%             while( (check1==0 ) && i<length(K) )
%                 if(DGfd.wd3(1,i) > tol ) 
%                     CFL_max = CFL(j-1);
%                     pp = strcat('DGp',num2str(P),' and RK',num2str(Prk),' Beta= ',num2str(Beta));
%                     disp(pp);
%                     pCFL = strcat('Unstable CFL: ',num2str(CFL(j))...
%                         ,', at K: ',num2str(K(i))...
%                         ,', |G|= ',num2str(exp(CFL(j)*DGfd.wd3(1,i))),', Omega1: ',num2str(DGfd.wd3(1,i)));
%                     disp(pCFL);
%                     maxCFL = strcat('CFL_max: ',num2str(CFL_max));
%                     disp(maxCFL);
%                     K_unstable_cfl=K(i);
%                     check1=1; 
%                 end
%                 i=i+1;
%             end 
%             
%             while( (check1==0 ) && i<length(K) )
%                 if(DGfd.wd4(1,i) > tol ) 
%                     CFL_max = CFL(j-1);
%                     pp = strcat('DGp',num2str(P),' and RK',num2str(Prk),' Beta= ',num2str(Beta));
%                     disp(pp);
%                     pCFL = strcat('Unstable CFL: ',num2str(CFL(j))...
%                         ,', at K: ',num2str(K(i))...
%                         ,', |G|= ',num2str(exp(CFL(j)*DGfd.wd4(1,i))),', Omega1: ',num2str(DGfd.wd4(1,i)));
%                     disp(pCFL);
%                     maxCFL = strcat('CFL_max: ',num2str(CFL_max));
%                     disp(maxCFL);
%                     K_unstable_cfl=K(i);
%                     check1=1; 
%                 end
%                 i=i+1;
%             end 
%             
%             while( (check1==0 ) && i<length(K) )
%                 if(DGfd.wd5(1,i) > tol ) 
%                     CFL_max = CFL(j-1);
%                     pp = strcat('DGp',num2str(P),' and RK',num2str(Prk),' Beta= ',num2str(Beta));
%                     disp(pp);
%                     pCFL = strcat('Unstable CFL: ',num2str(CFL(j))...
%                         ,', at K: ',num2str(K(i))...
%                         ,', |G|= ',num2str(exp(CFL(j)*DGfd.wd5(1,i))),', Omega1: ',num2str(DGfd.wd5(1,i)));
%                     disp(pCFL);
%                     maxCFL = strcat('CFL_max: ',num2str(CFL_max));
%                     disp(maxCFL);
%                     K_unstable_cfl=K(i);
%                     check1=1; 
%                 end
%                 i=i+1;
%             end
%             
%         end
    
    else
        break;
    end
end