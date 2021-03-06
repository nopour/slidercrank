function H_matrix = Hfunc(A_2,E2,I2,Rho2,dq1,dq2,dq3,dq4,dq5,g,l1,l2,m1,m3,q1,q2,q3,q4,q5)
%HFUNC
%    H_MATRIX = HFUNC(A_2,E2,I2,RHO2,DQ1,DQ2,DQ3,DQ4,DQ5,G,L1,L2,M1,M3,Q1,Q2,Q3,Q4,Q5)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    25-Apr-2022 11:56:53

t2 = cos(q1);
t3 = cos(q2);
t4 = sin(q2);
t5 = q1+q2;
t6 = dq1.^2;
t7 = dq2.^2;
t8 = l2.^2;
t9 = l2.^3;
t11 = q2.*2.0;
t16 = 1.0./l2;
t17 = -q2;
t10 = t8.^2;
t12 = cos(t11);
t13 = t3.^2;
t14 = sin(t11);
t15 = sin(t5);
t18 = q1+t17;
t21 = A_2.*Rho2.*g.*t3.*t9.*3.5e+1;
t19 = cos(t18);
t20 = sin(t18);
t22 = A_2.*Rho2.*l1.*t6.*t9.*t20.*3.5e+1;
H_matrix = [(l1.*(g.*m1.*t2.*6.0+dq2.*dq4.*m3.*t19.*1.2e+1+l2.*m3.*t7.*t15.*6.0+l2.*m3.*t7.*t20.*6.0+m3.*q4.*t7.*t15.*6.0+m3.*q4.*t7.*t20.*6.0-dq2.*dq4.*m3.*cos(t5).*1.2e+1+l1.*m3.*t6.*sin(q1.*2.0).*6.0+A_2.*Rho2.*g.*l2.*t2.*1.2e+1+A_2.*Rho2.*t7.*t8.*t20.*6.0+A_2.*Rho2.*dq2.*dq4.*l2.*t19.*1.2e+1+A_2.*Rho2.*dq2.*dq3.*t8.*t20.*2.0-A_2.*Rho2.*dq2.*dq5.*t8.*t20.*2.0+A_2.*Rho2.*l2.*q4.*t7.*t20.*6.0-A_2.*Rho2.*q3.*t7.*t8.*t19+A_2.*Rho2.*q5.*t7.*t8.*t19))./1.2e+1;dq2.*dq4.*l2.*m3+dq2.*dq4.*m3.*q4+(m3.*t7.*t8.*t14)./2.0+(m3.*q4.^2.*t7.*t14)./2.0+A_2.*Rho2.*dq2.*dq4.*t8.*(2.0./3.0)+(A_2.*Rho2.*g.*t3.*t8)./2.0-dq2.*dq4.*l2.*m3.*t12-dq2.*dq4.*m3.*q4.*t12+(l1.*l2.*m3.*t6.*t15)./2.0-(l1.*l2.*m3.*t6.*t20)./2.0+(l1.*m3.*q4.*t6.*t15)./2.0+l2.*m3.*q4.*t7.*t14-(l1.*m3.*q4.*t6.*t20)./2.0+A_2.*Rho2.*dq2.*dq4.*l2.*q4.*(2.0./3.0)+A_2.*Rho2.*dq2.*dq3.*q3.*t9.*(2.0./1.05e+2)-(A_2.*Rho2.*dq2.*dq3.*q5.*t9)./7.0e+1-(A_2.*Rho2.*dq2.*dq5.*q3.*t9)./7.0e+1+A_2.*Rho2.*dq2.*dq5.*q5.*t9.*(2.0./1.05e+2)+(A_2.*Rho2.*g.*l2.*q4.*t3)./2.0-(A_2.*Rho2.*g.*q3.*t4.*t8)./1.2e+1+(A_2.*Rho2.*g.*q5.*t4.*t8)./1.2e+1-(A_2.*Rho2.*l1.*t6.*t8.*t20)./2.0-(A_2.*Rho2.*l1.*l2.*q4.*t6.*t20)./2.0+(A_2.*Rho2.*l1.*q3.*t6.*t8.*t19)./1.2e+1-(A_2.*Rho2.*l1.*q5.*t6.*t8.*t19)./1.2e+1;(t16.*(t21-t22+E2.*I2.*q3.*1.68e+3+E2.*I2.*q5.*8.4e+2+A_2.*Rho2.*dq2.*dq4.*t9.*2.8e+1-A_2.*Rho2.*q3.*t7.*t10.*4.0+A_2.*Rho2.*q5.*t7.*t10.*3.0))./4.2e+2;t16.*(A_2.*E2.*q4.*-3.0e+1+m3.*t7.*t8.*t13.*3.0e+1+A_2.*Rho2.*t7.*t9.*1.0e+1+A_2.*Rho2.*dq2.*dq3.*t9.*2.0-A_2.*Rho2.*dq2.*dq5.*t9.*3.0-A_2.*Rho2.*g.*t4.*t8.*1.5e+1+A_2.*Rho2.*q4.*t7.*t8.*1.0e+1+dq2.*dq4.*l2.*m3.*t14.*3.0e+1+l2.*m3.*q4.*t7.*t13.*3.0e+1+l1.*l2.*m3.*t2.*t3.*t6.*3.0e+1+A_2.*Rho2.*l1.*t2.*t3.*t6.*t8.*1.5e+1+A_2.*Rho2.*l1.*t4.*t6.*t8.*sin(q1).*1.5e+1).*(-1.0./3.0e+1);(t16.*(-t21+t22+E2.*I2.*q3.*8.4e+2+E2.*I2.*q5.*1.68e+3-A_2.*Rho2.*dq2.*dq4.*t9.*4.2e+1+A_2.*Rho2.*q3.*t7.*t10.*3.0-A_2.*Rho2.*q5.*t7.*t10.*4.0))./4.2e+2];
