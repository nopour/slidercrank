clc

close all

addpath('GFunction');

global  l1 l2 m1 m2 m3 g Rho1 Rho2 A_1 A_2 E1 E2 I1 I2

% l1=6;
% l2=12;
% m1=1.6;
% m2=3.2;
% m3=0.0;
% g=32;


% l1=0.5;
% l2=1;
% m1=0.1832;
% m2=0.2;
% m3=0.1832;
% g=9.8;

l1=0.1524;
l2=0.304;
% m2=1;
m3=0.0;
g=0;

A_1=7.854E-5;
A_2=A_1;
Rho1=2.77E3;
Rho2=Rho1;
E1=1e9;
E2=0.5E8;
I1=4.909E-10;
I2=I1;
m1=A_1*Rho1*l1;
m2=A_2*Rho2*l2;


dt=0.01;
timeSpan=0:dt:1.6;
tsize=length(timeSpan);


q_0=[0;0];
% q_0=[0;0;0;0;0;0;0;0];
dq_0=[0;0];
mu_0=0;

options=odeset('maxstep',1e-4);
%%
t0=clock;
[~,zIMM]=ode45(@IMM_dynamics,timeSpan,[q_0;dq_0;mu_0],options);
t1=clock;
timeIMM=etime(t1,t0);
disp(['IMM sim time: ' num2str(timeIMM) '(s)'])

%%
t0=clock;
[~,zAM]=ode45(@AM_dynamics,timeSpan,[q_0;dq_0],options);
t1=clock;
timeAM=etime(t1,t0);
disp(['AM sim time: ' num2str(timeAM) '(s)'])

%%
t0=clock;
[~,zEM]=ode45(@EM_dynamics,timeSpan,[q_0;dq_0],options);
t1=clock;
timeEM=etime(t1,t0);
disp(['EM sim time: ' num2str(timeEM) '(s)'])

%%
t0=clock;
[~,zGM]=ode45(@GM_dynamics,timeSpan,[q_0;dq_0],options);
t1=clock;
timeGM=etime(t1,t0);
disp(['GM sim time: ' num2str(timeGM) '(s)'])

%%
IM.q1=zIMM(:,1);
IM.q2=zIMM(:,2);

IM.dq1=zIMM(:,3);
IM.dq2=zIMM(:,4);

%
AM.q1=zAM(:,1);
AM.q2=zAM(:,2);

AM.dq1=zAM(:,3);
AM.dq2=zAM(:,4);

%
EM.q1=zEM(:,1);
EM.q2=zEM(:,2);

EM.dq1=zEM(:,3);
EM.dq2=zEM(:,4);

%
GM.q1=zGM(:,1);
GM.q2=zGM(:,2);

GM.dq1=zGM(:,3);
GM.dq2=zGM(:,4);

%% 
%constraint
AE_IM = AEfunc(IM.dq1,IM.dq2,l1,l2,IM.q1,IM.q2);
AE_AM = AEfunc(AM.dq1,AM.dq2,l1,l2,AM.q1,AM.q2);
AE_EM = AEfunc(EM.dq1,EM.dq2,l1,l2,EM.q1,EM.q2);
AE_GM = AEfunc(GM.dq1,GM.dq2,l1,l2,GM.q1,GM.q2);

% Energy

E_IM=Efunc(IM.dq1,IM.dq2,g,l1,l2,m1,m2,m3,IM.q1,IM.q2);
E_AM=Efunc(AM.dq1,AM.dq2,g,l1,l2,m1,m2,m3,AM.q1,AM.q2);
E_EM=Efunc(EM.dq1,EM.dq2,g,l1,l2,m1,m2,m3,EM.q1,EM.q2);
E_GM=Efunc(GM.dq1,GM.dq2,g,l1,l2,m1,m2,m3,GM.q1,GM.q2);


ER_IM=(E_IM-E_IM(1))/E_IM(1);
ER_AM=(E_IM-E_AM(1))/E_AM(1);
ER_EM=(E_EM-E_EM(1))/E_EM(1);
ER_GM=(E_GM-E_GM(1))/E_GM(1);

%%
%figure
%hold on; grid on
%plot(timeSpan,ER_IM,'r-')
%plot(timeSpan,ER_AM,'g-')
%plot(timeSpan,ER_EM,'b-')
%plot(timeSpan,ER_GM,'k-')
%legend('IM','AM','EM','GM')
%xlabel('time(s)')
%title('Mechanical Energy Error')

figure
hold on; grid on
plot(timeSpan,AE_IM,'r-','linewidth',1.5)
plot(timeSpan,AE_AM,'g.','linewidth',2.5)
plot(timeSpan,AE_EM,'b-.','linewidth',1.5)
plot(timeSpan,AE_GM,'k--','linewidth',1.5)
legend('IM','AM','EM','GM')
xlabel('time(s)')
grid minor
title('Constraint Error')
%%
figure
hold on; grid on
plot(timeSpan,(IM.q1)*180/pi,'r-','linewidth',2.5)
plot(timeSpan,(AM.q1)*180/pi,'g--','linewidth',2)
plot(timeSpan,(EM.q1)*180/pi,'b-','linewidth',1.5)
plot(timeSpan,(GM.q1)*180/pi,'k-','linewidth',1)
legend('IM','AM','EM','GM')
xlabel('time(s)')
title('\theta_1')


figure
hold on; grid on
plot(timeSpan,(IM.q2)*180/pi,'r-','linewidth',2.5)
plot(timeSpan,(AM.q2)*180/pi,'g--','linewidth',2)
plot(timeSpan,(EM.q2)*180/pi,'b-','linewidth',1.5)
plot(timeSpan,(GM.q2)*180/pi,'k-','linewidth',1)
legend('IM','AM','EM','GM')
xlabel('time(s)')
title('\theta_2')
ylabel('Angle (deg)')

%%

figure
% A=openfig('untitled.fig');
hold on
sl=(l1)*cos(GM.q1)+(l2)*cos(GM.q2);
plot(timeSpan,sl,'-k','LineWidth',1.5);
% legend('Flex','Rigid');
grid minor
xlabel('time(s)')
ylabel('slider position(m)')
ylim([0.1 0.5])

%%
% %
% %%
% 
% figure 
% for i=1:length(timeSpan)
%     
%     XO(i)=0;
%     YO(i)=0;
%     
%     XA(i)=XO(i)+l1*cos(AM.q1(i));
%     YA(i)=YO(i)+l1*sin(AM.q1(i));
%     
%     XB(i)=XA(i)+l2*cos(AM.q2(i));
%     YB(i)=YA(i)+l2*sin(AM.q2(i));
%     
%     P1x=[XO(i);XA(i)];
%     P1y=[YO(i);YA(i)];
%     
%     P2x=[XA(i);XB(i)];
%     P2y=[YA(i);YB(i)];
%     
%     
% %     plot(xR,yR,'y-','linewidth',4)
% %     hold on
% %     plot(XO(i)+xd,YO(i)+yd,'r-','linewidth',4)
%     
%     plot(P1x,P1y,'b-','linewidth',1)
%     hold on
%     plot(P2x,P2y,'g-','linewidth',2)
%     
%     axis([-l1-l2 l1+l2 -(l1+l2) l1+l2])
%     grid minor
% 
%     axis equal
%     
%     hold off
%     pause(0.01)
%     
% end
    
