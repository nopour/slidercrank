function AE = AEfunc(dq1,dq2,l1,l2,q1,q2)
%AEFUNC
%    AE = AEFUNC(DQ1,DQ2,L1,L2,Q1,Q2)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    28-Apr-2022 07:39:39

AE = (dq1.*l1.*cos(q1)+dq2.*l2.*cos(q2)).^2;
