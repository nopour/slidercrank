function dZ=EM_dynamics(t,Z)

dZ=Z;

global A_1 A_2 Rho1 Rho2 l1 l2 m3 E1 E2 I1 I2 g
 

q=Z(1:8); dq=Z(9:16);

q1=q(1); q2=q(2); q3=q(3); q4=q(4); 
q5=q(5); q6=q(6); q7=q(7); q8=q(8); 

dq1=dq(1); 
dq2=dq(2); 
dq3=dq(3); 
dq4=dq(4); 
dq5=dq(5); 
dq6=dq(6); 
dq7=dq(7); 
dq8=dq(8); 

M = Mfunc(A_1,A_2,Rho1,Rho2,l1,l2,m3,q1,q2,q3,q4,q5,q6,q7,q8);
H = Hfunc(A_1,A_2,E1,E2,I1,I2,Rho1,Rho2,dq1,dq2,dq3,dq4,dq5,dq6,dq7,dq8,g,l1,l2,m3,q1,q2,q3,q4,q5,q6,q7,q8);
A1 =A1func(l1,l2,q1,q2,q4,q7).';
A1dot = A1dotfunc(dq1,dq2,dq4,dq7,l1,l2,q1,q2,q4,q7).';

if t<=0.7
    TO=0.01*(1-exp(-t/0.167));
elseif t>=0.7
    TO=0;
end

FTT=[TO;0;TO;0;0;0;0;0];


F=FTT-H;

M1=M(1:1,:);  M2=M(2:8,:);
F1=F(1:1,:);  F2=F(2:8,:);
a1=A1(:,1:1);  a2=A1(:,2:8);


M_E=[(M2-a2'*(pinv(a1')*M1));-A1];

F_E=[(F2-a2'*(pinv(a1')*F1));A1dot*dq];

ddq=pinv(M_E)*F_E;

dZ=[dq;ddq];

end