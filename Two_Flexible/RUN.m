clc
clear
close all

addpath('GFunction');

global  l1 l2 m3 g Rho1 Rho2 A_1 A_2 E1 E2 I1 I2
% 
l1=0.1524;
l2=0.304;
% m2=1;

g=0;

A_1=7.854E-5;
A_2=A_1;
Rho1=2.77E3;
Rho2=Rho1;
E1=1e9;
E2=0.5E8;
I1=4.909E-10;
I2=I1;
m3=0;

% 
% l1=6;
% l2=12;
% m3=0;
% g=0;
% A_1=0.25^2/4*pi;
% A_2=A_1;
% Rho1=0.0007331;
% Rho2=Rho1;      %similar %$
% E1=30*10^6;
% E2=E1;
% I1=0.0001917476;
% I2=I1;
% % m1=A_1*l1*Rho1;




dt=0.0001;
timeSpan=0:dt:1.6;
tsize=length(timeSpan);



q_0=[0;0;0;0;0;0;0;0];
dq_0=[0;0;0;0;0;0;0;0];
mu_0=[0];

options=odeset('maxstep',1e-3);
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
IM.q3=zIMM(:,3);
IM.q4=zIMM(:,4);
IM.q5=zIMM(:,5);
IM.q6=zIMM(:,6);
IM.q7=zIMM(:,7);
IM.q8=zIMM(:,8);

IM.dq1=zIMM(:,9);
IM.dq2=zIMM(:,10);
IM.dq3=zIMM(:,11);
IM.dq4=zIMM(:,12);
IM.dq5=zIMM(:,13);
IM.dq6=zIMM(:,14);
IM.dq7=zIMM(:,15);
IM.dq8=zIMM(:,16);

%
AM.q1=zAM(:,1);
AM.q2=zAM(:,2);
AM.q3=zAM(:,3);
AM.q4=zAM(:,4);
AM.q5=zAM(:,5);
AM.q6=zAM(:,6);
AM.q7=zAM(:,7);
AM.q8=zAM(:,8);

AM.dq1=zAM(:,9);
AM.dq2=zAM(:,10);
AM.dq3=zAM(:,11);
AM.dq4=zAM(:,12);
AM.dq5=zAM(:,13);
AM.dq6=zAM(:,14);
AM.dq7=zAM(:,15);
AM.dq8=zAM(:,16);
%
EM.q1=zEM(:,1);
EM.q2=zEM(:,2);
EM.q3=zEM(:,3);
EM.q4=zEM(:,4);
EM.q5=zEM(:,5);
EM.q6=zEM(:,6);
EM.q7=zEM(:,7);
EM.q8=zEM(:,8);

EM.dq1=zEM(:,9);
EM.dq2=zEM(:,10);
EM.dq3=zEM(:,11);
EM.dq4=zEM(:,12);
EM.dq5=zEM(:,13);
EM.dq6=zEM(:,14);
EM.dq7=zEM(:,15);
EM.dq8=zEM(:,16);
%
GM.q1=zGM(:,1);
GM.q2=zGM(:,2);
GM.q3=zGM(:,3);
GM.q4=zGM(:,4);
GM.q5=zGM(:,5);
GM.q6=zGM(:,6);
GM.q7=zGM(:,7);
GM.q8=zGM(:,8);

GM.dq1=zGM(:,9);
GM.dq2=zGM(:,10);
GM.dq3=zGM(:,11);
GM.dq4=zGM(:,12);
GM.dq5=zGM(:,13);
GM.dq6=zGM(:,14);
GM.dq7=zGM(:,15);
GM.dq8=zGM(:,16);
%% 
%constraint
AE_IM = AEfunc(IM.dq1,IM.dq2,IM.dq4,IM.dq7,l1,l2,IM.q1,IM.q2,IM.q4,IM.q7);
AE_AM = AEfunc(AM.dq1,AM.dq2,AM.dq4,AM.dq7,l1,l2,AM.q1,AM.q2,AM.q4,AM.q7);
AE_EM = AEfunc(EM.dq1,EM.dq2,EM.dq4,EM.dq7,l1,l2,EM.q1,EM.q2,EM.q4,EM.q7);
AE_GM = AEfunc(GM.dq1,GM.dq2,GM.dq4,GM.dq7,l1,l2,GM.q1,GM.q2,GM.q4,GM.q7);

%%
figure
hold on; grid on
plot(timeSpan,AE_IM,'r-')
plot(timeSpan,AE_AM,'g-')
plot(timeSpan,AE_EM,'b-')
plot(timeSpan,AE_GM,'k-')
legend('IM','AM','EM','GM')
xlabel('time(s)')
title('Constraint Error')
%%
% figure
% hold on; grid on
% plot(timeSpan,(IM.q1)*180/pi,'r-','linewidth',2.5)
% plot(timeSpan,(AM.q1)*180/pi,'g--','linewidth',2)
% plot(timeSpan,(EM.q1)*180/pi,'b-','linewidth',1.5)
% plot(timeSpan,(GM.q1)*180/pi,'k-','linewidth',1)
% legend('IM','AM','EM','GM')
% xlabel('time(s)')
% title('\theta_1')

% 
% figure
% hold on; grid on
% plot(timeSpan,(IM.q2)*180/pi,'r-','linewidth',2.5)
% plot(timeSpan,(AM.q2)*180/pi,'g--','linewidth',2)
% plot(timeSpan,(EM.q2)*180/pi,'b-','linewidth',1.5)
% plot(timeSpan,(GM.q2)*180/pi,'k-','linewidth',1)
% legend('IM','AM','EM','GM')
% xlabel('time(s)')
% title('\theta_2')
% 
% 
% figure
% hold on; grid on
% plot(timeSpan,IM.q4,'r-','linewidth',2.5)
% plot(timeSpan,AM.q4,'g--','linewidth',2)
% plot(timeSpan,EM.q4,'b-','linewidth',1.5)
% plot(timeSpan,GM.q4,'k-','linewidth',1)
% legend('IM','AM','EM','GM')
% xlabel('time(s)')
% title('\delta_2')

%%
% figure
% hold on; grid on
% plot(timeSpan,IM.q6,'r-','linewidth',2.5)
% plot(timeSpan,AM.q6,'g--','linewidth',2)
% plot(timeSpan,EM.q6,'b-','linewidth',1.5)
% plot(timeSpan,GM.q6,'k-','linewidth',1)
% legend('IM','AM','EM','GM')
% xlabel('time(s)')
% title('\delta_2')

%%
% SL1=(l1+AM.q4).*sin(AM.q1)+(l2+AM.q7).*sin(AM.q2);
SL2=(l1+IM.q4).*cos(IM.q1)+(l2+IM.q7).*cos(IM.q2);

figure 
% plot(timeSpan,SL1,'r','linewidth',1.5)
hold on
plot(timeSpan,SL2,'k--','linewidth',2)
axis([0 1.6 0.1 0.5])
grid minor

%%
% figure 
% for i=1:length(timeSpan)
%     
%     XO(i)=0;
%     YO(i)=0;
%     
%     XA(i)=XO(i)+(l1+IM.q4(i))*cos(AM.q1(i));
%     YA(i)=YO(i)+(l1+IM.q4(i))*sin(AM.q1(i));
%     
%     XB(i)=XA(i)+(l2+IM.q7(i))*cos(AM.q2(i));
%     YB(i)=YA(i)+(l2+IM.q7(i))*sin(AM.q2(i));
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
%     pause(0.00001)
%     
% end
% %     
