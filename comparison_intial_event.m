clc;clear;close all;
%tarif soorat soal
A1=5;A2=4;A3=5;A4=6;E=36000;sigmay=25;alpha=0.03;    A=[5,4,5,6];A_n=diag(A);

L1=sqrt((9^2)+(5^2));L2=sqrt((9^2)+(17^2));L3=sqrt((8^2)+(17^2));L4=sqrt((15^2)+(6^2));L=[L1,L2,L3,L4];

teta1=3.65;teta2=4.225;teta3=5.152;teta4=5.903;%rad

qy=[sigmay*A1,sigmay*A2,sigmay*A3,sigmay*A4];

a=[-cos(teta1),-sin(teta1);-cos(teta2),-sin(teta2);-cos(teta3),-sin(teta3);-cos(teta4),-sin(teta4)];

Kel=zeros(4);                Kel(1,1)=A1*E/L1;        Kel(2,2)=A2*E/L2;       Kel(3,3)=A3*E/L3 ;      Kel(4,4)=A4*E/L4;

a_s=zeros(4);              a_s(1,1)=1/L1;        a_s(2,2)=1/L2;       a_s(3,3)=1/L3 ;      a_s(4,4)=1/L4;

Disp=[0;0];

p_Ext=[250;208.3];

tol=1e-10;

%tavabe komaki tarsim
x1=0:0.001:0.035;
y1=p_Ext(1,1).*(x1>=0);
y2=p_Ext(2,1).*(x1>=0);
y=sqrt((p_Ext(1,1)^2)+(p_Ext(2,1)^2)).*(x1>=0);
n=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% method full newton 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
step=input('enter how many steps do you want to use in calculating:');
p_ext=p_Ext./step
for m1=1:step
    m2=1;
         if m1==1
              p_r=m1*p_ext;
         else 
             p_r=m1*p_ext-(m1-1)*p_ext;
         end

          

    
            while m2<=100 & p_r(1,1)>=tol &p_r(2,1)>=tol
                     
                     
                p_int=[0;0];
                  Delta_u=inv(a'*Kel*a)*p_r;
                  disp_u=Disp+Delta_u;
                  Disp=disp_u;
                  V_n=a*disp_u;
       
                  epsilon=a_s*V_n;
                  sigma_n=zeros(4,1);
                  for j=1:4
                      if epsilon(j,1)>=6.95e-4;
                         sigma_n(j,1)=(epsilon(j,1)-6.95e-4)*E*alpha+25;
                      else 
                         sigma_n(j,1)=epsilon(j,1)*E;
                      end
                  end
                  q_n=A_n*sigma_n;
                  p_int=a'*q_n;
                    p_r=m1*p_ext-p_int;
                  
                    n=n+1;
                    m2=m2+1;
                    
                      
                    
                  
                  b_n(1:2,n)=p_int;
                  c_n(1:2,n)=disp_u;
                  
                  subplot(3,1,1)
                  title('p1-u1')
                  plot(c_n(1,n),b_n(1,n),'ob',x1,y1),grid on,hold on
                  b_nf=b_n(1,1:n);
                  c_nf=c_n(1,1:n);
                  line(c_nf,b_nf)
                 % 
                   subplot(3,1,2)
                  title('p2-u2')
                  plot(c_n(2,n),b_n(2,n),'ob',x1,y2),grid on,hold on
                  b_nfu=b_n(2,1:n);
                  c_nfu=c_n(2,1:n);
                  line(c_nfu,b_nfu)
                  
                  B_n(1,n)=sqrt(((b_n(1,n))^2)+((b_n(2,n))^2));
                  C_n(1,n)=sqrt(((c_n(1,n))^2)+((c_n(2,n))^2));
                   subplot(3,1,3)
                  title('p-u')
                  plot(C_n,B_n,'ob',x1,y),grid on,hold on
                  B_nf=B_n(1,1:n);
                  C_nf=C_n(1,1:n);
                  line(C_nf,B_nf)
                 
                     end
               
       
end 

                  
                   
                   
                  
%tarif soorat soal
A1=5;A2=4;A3=5;A4=6;E=36000;sigmay=25;
A=[5,4,5,6];
L1=sqrt((9^2)+(5^2));L2=sqrt((9^2)+(17^2));L3=sqrt((8^2)+(17^2));L4=sqrt((15^2)+(6^2));
L=[L1,L2,L3,L4];
teta1=3.65;teta2=4.225;teta3=5.152;teta4=5.903;%rad
qy=[sigmay*A1,sigmay*A2,sigmay*A3,sigmay*A4];

a=[-cos(teta1),-sin(teta1);-cos(teta2),-sin(teta2);-cos(teta3),-sin(teta3);-cos(teta4),-sin(teta4)];
Kel=zeros(4);                Kel(1,1)=A1*E/L1;        Kel(2,2)=A2*E/L2;       Kel(3,3)=A3*E/L3 ;      Kel(4,4)=A4*E/L4;
    
    landa=1;
    p=landa*[3,2.5] ;
     p_event=[0;0];
  U_event=[0;0];
  q_event=[0;0;0;0];
  landaprim=0;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 %event_to_event
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
    b=zeros(2,5);
    c=zeros(2,5);
 
    for i=1:5         %shomarande event
            
          landa=1;
          p=landa*p;
          U=inv(a'*Kel*a)*p';
          V=a*U;
          q=Kel*V;
          if i==1
              Delta_q=qy;
          else    
          Delta_q=qy'-abs(q_event);
          end
                 for k=1:4
                      if Delta_q(k)<=0 
                      landaprime(k)=1000;
                      else
                     landaprime(k)=Delta_q(k)/q(k);
                      end
                 end
                     landaprim=min(abs(landaprime));
                  
                        Q=landaprim*q ;
                        u=landaprim*U;
                        P=landaprim*p;
                        
                         b(1:2,i)=p_event;
                         c(1:2,i)=U_event ;               
                        
                        p_event=P'+p_event;
                         U_event=u+U_event;
                         q_event=Q+q_event;
                         %
                        
                       
               for j=1:4
                   if abs(q_event(j))==qy(j);
                     Kel(j,j)=0.03*Kel(j,j);
                   end
               end
    end
    
    w=zeros(1,5);
    z=zeros(1,5);
    for j=1:5
   
    subplot(3,1,1)
    title('p1-U1')
   
    plot(c(1,j),b(1,j),'*r');hold on;grid on
    x=b(1,1:j);
    y=c(1,1:j);
    line(y,x,'color','r')
    
   subplot(3,1,2)
    title('p2-U2')
    plot(c(2,j),b(2,j),'*r');hold on;grid on
    e=b(2,1:j);
    f=c(2,1:j);
    line(f,e,'color','r')
    
     subplot(3,1,3)
    title('p-U')
    B=sqrt(((b(1,j))^2)+((b(2,j))^2));
    C=sqrt(((c(1,j))^2)+((c(2,j))^2));
    w(1,j)=B;
    z(1,j)=C;
    plot(z(1,j),w(1,j),'*r');hold on;grid on
    g=w(1,1:j);
    h=z(1,1:j);
    line(h,g,'color','r')
    end
        
    


          
          
         
    
    

    
    
    