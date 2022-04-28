function A1dot = A1dotfunc(dq1,dq2,dq4,dq7,l1,l2,q1,q2,q4,q7)
%A1DOTFUNC
%    A1DOT = A1DOTFUNC(DQ1,DQ2,DQ4,DQ7,L1,L2,Q1,Q2,Q4,Q7)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    28-Apr-2022 15:26:39

t2 = cos(q1);
t3 = cos(q2);
A1dot = [dq4.*t2-dq1.*sin(q1).*(l1+q4);dq7.*t3-dq2.*sin(q2).*(l2+q7);0.0;dq1.*t2;0.0;0.0;dq2.*t3;0.0];