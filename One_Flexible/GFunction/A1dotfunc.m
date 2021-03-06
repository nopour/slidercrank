function A1dot = A1dotfunc(dq1,dq2,dq4,l1,l2,q1,q2,q4)
%A1DOTFUNC
%    A1DOT = A1DOTFUNC(DQ1,DQ2,DQ4,L1,L2,Q1,Q2,Q4)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    25-Apr-2022 11:56:54

t2 = cos(q2);
A1dot = [-dq1.*l1.*sin(q1);dq4.*t2-dq2.*sin(q2).*(l2+q4);0.0;dq2.*t2;0.0];
