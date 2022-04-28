function AE = AEfunc(dq1,dq2,dq4,dq7,l1,l2,q1,q2,q4,q7)
%AEFUNC
%    AE = AEFUNC(DQ1,DQ2,DQ4,DQ7,L1,L2,Q1,Q2,Q4,Q7)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    28-Apr-2022 15:26:44

AE = (dq4.*sin(q1)+dq7.*sin(q2)+dq1.*cos(q1).*(l1+q4)+dq2.*cos(q2).*(l2+q7)).^2;
