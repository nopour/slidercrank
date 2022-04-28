function dZ=GM_dynamics(t,Z)

dZ=Z;


global m1 m2 l1 l2 m3 g
 

q=Z(1:2); dq=Z(3:4);

q1=q(1); q2=q(2); 

dq1=dq(1); 
dq2=dq(2); 



M = Mfunc(l1,l2,m1,m2,m3,q1,q2);
H = Hfunc(dq1,dq2,g,l1,l2,m1,m2,m3,q1,q2);
A1 = A1func(l1,l2,q1,q2).';
A1dot = A1dotfunc(dq1,dq2,l1,l2,q1,q2).';



% if t<=0.7
    TO=0.01*(1-exp(-t/0.167));
% elseif t>=0.7 
%     TO=0;
% end

FTT=[TO;0];
F=FTT-H;


F_G=F-A1'*(pinv(A1*(pinv(M)*A1'))*(A1dot*dq+A1*(pinv(M)*F)));

ddq=pinv(M)*F_G;

dZ=[dq;ddq];


end