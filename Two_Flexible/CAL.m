clc
clear all
close all
%%
% q1=theta1   q2=theta2 

addpath('GFunction')

syms t x
syms q1 q2 q3 q4 q5 q6 q7 q8
syms dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8
syms l1 l2 m1 m2 m3 g Rho1 A_1 E1 I1 Rho2 A_2 E2 I2;
tor1=0;


%% Link 1
Q1=[cos(q1) -sin(q1);sin(q1) cos(q1)];
z1=x/l1;
S1=[0 z1 0;l1*(z1^3-2*z1^2+z1) 0 l1*(z1^3-z1^2)];
t1=[q3 q4 q5].';
r1(x)=Q1*([x;0]+S1*t1);
v1(x)=diff(r1,q1)*dq1+diff(r1,q3)*dq3+diff(r1,q4)*dq4+diff(r1,q5)*dq5;
T1=0.5*Rho1*A_1*int(v1.'*v1,x,0,l1);
u1=S1*t1;
VE1=0.5*E1*I1*int((diff(u1.'*[0;1],x,2))^2,x,0,l1) ...
    +0.5*E1*A_1*int((diff(u1.'*[1;0],x,1))^2,x,0,l1);
VG1=Rho1*A_1*g*int(r1.'*[0;1],x,0,l1);
V1=VE1+VG1;
%% Link 2
Q2=[cos(q2) -sin(q2);sin(q2) cos(q2)];
z2=x/l2;
S2=[0 z2 0;l2*(z2^3-2*z2^2+z2) 0 l2*(z2^3-z2^2)];
t2=[q6 q7 q8].';
r2(x)=r1(l1)+Q2*([x;0]+S2*t2);
v2(x)=diff(r2,q1)*dq1+diff(r2,q2)*dq2+diff(r2,q3)*dq3+diff(r2,q4)*dq4+diff(r2,q5)*dq5 ...
     +diff(r2,q6)*dq6+diff(r2,q7)*dq7+diff(r2,q8)*dq8;
T2=0.5*Rho2*A_2*int(v2.'*v2,x,0,l2);
u2=S2*t2;
VE2=0.5*E2*I2*int((diff(u2.'*[0;1],x,2))^2,x,0,l2) ...
    +0.5*E2*A_2*int((diff(u2.'*[1;0],x,1))^2,x,0,l2);
VG2=Rho2*A_2*g*int(r2.'*[0;1],x,0,l2);
V2=VE2+VG2;

%% Link 3
v3=v2(l2).'*[1;0];
T3=0.5*m3*v3^2;
%% Energy
T=T1+T2+T3;
V=V1+V2;
E=T+V;

q=[q1 q2 q3 q4 q5 q6 q7 q8].';
dq=[dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8].';

Mass_matrix=jacobian(jacobian(T,dq),dq);
H_matrix=(jacobian(jacobian(T,dq),q))*dq-jacobian(T,q).'+jacobian(V,q).';
M1=[tor1;0;tor1;0;0;0;0;0];

A=(l1+q4)*sin(q1)+(l2+q7)*sin(q2);
A=diff(A,q1)*dq1+diff(A,q2)*dq2+diff(A,q3)*dq3+diff(A,q4)*dq4+diff(A,q5)*dq5 ...
     +diff(A,q6)*dq6+diff(A,q7)*dq7+diff(A,q8)*dq8;
 
A1=jacobian(A,dq).';
AE=simplify(A.'*A);


A1dot=A1;
for n=1:size(A1,1)
       A1dot(n,:)= (jacobian(A1(n,:),q)*dq).';
end

M=simplify(Mass_matrix);
H_matrix=simplify(H_matrix);
A1=simplify(A1);
A1dot=simplify(A1dot);
%%

Mfunc=matlabFunction(M,'File','GFunction\Mfunc');
Hfunc=matlabFunction(H_matrix,'File','GFunction\Hfunc');
A1func=matlabFunction(A1,'File','GFunction\A1func');
A1dotfunc=matlabFunction(A1dot,'File','GFunction\A1dotfunc');
Efunc=matlabFunction(simplify(E),'File','GFunction\Efunc');
AEfunc=matlabFunction((AE),'File','GFunction\AEfunc');

