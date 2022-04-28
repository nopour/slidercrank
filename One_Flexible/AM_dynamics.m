function dZ=AM_dynamics(t,Z)

dZ=Z;

global  A_2 Rho2 l1 l2 m3 m1 E2  I2 g
 

q=Z(1:5); dq=Z(6:10);

q1=q(1); q2=q(2); q3=q(3); q4=q(4); 
q5=q(5); 

dq1=dq(1); 
dq2=dq(2); 
dq3=dq(3); 
dq4=dq(4); 
dq5=dq(5); 


M = Mfunc(A_2,Rho2,l1,l2,m1,m3,q1,q2,q3,q4,q5);
H = Hfunc(A_2,E2,I2,Rho2,dq1,dq2,dq3,dq4,dq5,g,l1,l2,m1,m3,q1,q2,q3,q4,q5);
A1 = A1func(l1,l2,q1,q2,q4).';
A1dot = A1dotfunc(dq1,dq2,dq4,l1,l2,q1,q2,q4).';

if t<=0.7
    TO=0.01*(1-exp(-t/0.167));
elseif t>=0.7
    TO=0;
end


FTT=[TO;0;0;0;0];

F=FTT-H;

M_AM=[M -A1.';-A1 zeros(1)];
F_AM=[F;A1dot*dq];

X=pinv(M_AM)*F_AM;

dZ=[dq;X(1:5)];
end