function A1 = A1func(l1,l2,q1,q2,q4,q7)
%A1FUNC
%    A1 = A1FUNC(L1,L2,Q1,Q2,Q4,Q7)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    28-Apr-2022 15:26:39

A1 = [cos(q1).*(l1+q4);cos(q2).*(l2+q7);0.0;sin(q1);0.0;0.0;sin(q2);0.0];
