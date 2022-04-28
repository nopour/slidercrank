function dZ=EM_dynamics(t,Z)

dZ=Z;

global A_1 A_2 Rho1 Rho2 l1 l2 m3 E1 E2 I1 I2 g m1
 

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

M1=M(1:1,:);  M2=M(2:5,:);
F1=F(1:1,:);  F2=F(2:5,:);
a1=A1(:,1:1);  a2=A1(:,2:5);


M_E=[(M2-a2'*(pinv(a1')*M1));-A1];

F_E=[(F2-a2'*(pinv(a1')*F1));A1dot*dq];

ddq=pinv(M_E)*F_E;

dZ=[dq;ddq];

end