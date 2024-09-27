clc;clear;close all
A=[5,4,5,6];A1=5;A2=4;A3=5;A4=6;E=36000;sigmay=25;alpha=0.03;    A=[5,4,5,6];A_n=diag(A);
L1=sqrt((9^2)+(5^2));L2=sqrt((9^2)+(17^2));L3=sqrt((8^2)+(17^2));L4=sqrt((15^2)+(6^2));
L=[L1,L2,L3,L4];
teta1=3.65;teta2=4.225;teta3=5.152;teta4=5.903;%rad
qy=[sigmay*A1,sigmay*A2,sigmay*A3,sigmay*A4];

a=[-cos(teta1),-sin(teta1);-cos(teta2),-sin(teta2);-cos(teta3),-sin(teta3);-cos(teta4),-sin(teta4)];
Kel=zeros(4);                Kel(1,1)=A1*E/L1;        Kel(2,2)=A2*E/L2;       Kel(3,3)=A3*E/L3 ;      Kel(4,4)=A4*E/L4;
a_s=zeros(4);              a_s(1,1)=1/L1;        a_s(2,2)=1/L2;       a_s(3,3)=1/L3 ;      a_s(4,4)=1/L4;





R=[250;208.3]
lamba=0
pint=[0;0]
pr=R-pint
delta_u_bar=inv(a'*Kel*a)*(lamba*R-pint)
delta_u_2bar=inv(a'*Kel*a)*R
u0=[0;0]



r0=0.35
beta=50


syms x delta_lamba
eq1=((lamba+x)^2)+(1/beta)*((u0+delta_u_bar+delta_u_2bar*x)'*(u0+delta_u_bar+delta_u_2bar*x))==(r0^2)
delta_lamba=solve(eq1,x)

lamba=lamba+delta_lamba

u=u0+delta_u_bar+delta_u_2bar*delta_lamba

v=a*u

epsilon=a_s*v

sigma=zeros(4,1);
                
                  
                  
                  
                  
                  for j=1:4
                      if epsilon(j,1)>=6.95e-4;
                         sigma(j,1)=(epsilon(j,1)-6.95e-4)*E*alpha+25;
                      else 
                         sigma(j,1)=epsilon(j,1)*E;
                      end
                  end
                  
                  
   q=sigma*A_n
   
   pint=a'*q
   
   pr=lamba*R-pint
   
   
                  
                  
                  
                  
                  
