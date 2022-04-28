function H_matrix = Hfunc(dq1,dq2,g,l1,l2,m1,m2,m3,q1,q2)
%HFUNC
%    H_MATRIX = HFUNC(DQ1,DQ2,G,L1,L2,M1,M2,M3,Q1,Q2)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    28-Apr-2022 07:39:38

t2 = cos(q1);
t3 = cos(q2);
t4 = sin(q2);
t5 = dq1.^2;
t6 = dq2.^2;
t7 = -q2;
t8 = q1+t7;
t9 = sin(t8);
H_matrix = [(l1.*(g.*m1.*t2+g.*m2.*t2.*2.0+l2.*m2.*t6.*t9+l2.*m3.*t6.*t9+l2.*m3.*t6.*sin(q1+q2)+l1.*m3.*t5.*sin(q1.*2.0)))./2.0;(g.*l2.*m2.*t3)./2.0+l2.^2.*m3.*t3.*t4.*t6+(l1.*l2.*m2.*t2.*t4.*t5)./2.0+l1.*l2.*m3.*t2.*t4.*t5-(l1.*l2.*m2.*t3.*t5.*sin(q1))./2.0];
