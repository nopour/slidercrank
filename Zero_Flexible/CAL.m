clc
clear all
close all
%%
% q1=theta1   q2=theta2 

addpath('GFunction')

syms t x
syms q1 q2 
syms dq1 dq2 
syms l1 l2 m1 m2 m3 g 
tor1=0;


%% Link 1
%% Link 1
Q1=[cos(q1) -sin(q1);sin(q1) cos(q1)];
r1(x)=Q1*[x;0];
v1(x)=diff(r1,q1)*dq1;
r1(l1);
v1(l1);
T1=0.5*m1/l1*int(v1.'*v1,x,0,l1);
V1=m1*g/l1*int(r1.'*[0;1],x,0,l1);

%% Link 2
Q2=[cos(q2) -sin(q2);sin(q2) cos(q2)];
r2(x)=r1(l1)+Q2*([x;0]);
v2(x)=diff(r2,q1)*dq1+diff(r2,q2)*dq2;
T2=0.5*m2/l2*int(v2.'*v2,x,0,l2);
V2=m2*g/l2*int(r2.'*[0;1],x,0,l2);

%% Link 3
v3=v2(l2).'*[1;0];
T3=0.5*m3*v3^2;

%% Energy
T=T1+T2+T3;
V=V1+V2;
E=T+V;

q=[q1 q2].';
dq=[dq1 dq2].';

Mass_matrix=jacobian(jacobian(T,dq),dq);
H_matrix=(jacobian(jacobian(T,dq),q))*dq-jacobian(T,q).'+jacobian(V,q).';
M1=[tor1;0];

A=(l1)*sin(q1)+(l2)*sin(q2);
A=diff(A,q1)*dq1+diff(A,q2)*dq2;
 
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

