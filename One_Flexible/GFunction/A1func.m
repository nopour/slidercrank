function A1 = A1func(l1,l2,q1,q2,q4)
%A1FUNC
%    A1 = A1FUNC(L1,L2,Q1,Q2,Q4)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    25-Apr-2022 11:56:53

A1 = [l1.*cos(q1);cos(q2).*(l2+q4);0.0;sin(q2);0.0];
