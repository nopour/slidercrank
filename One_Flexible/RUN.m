clc

close all

addpath('GFunction');

global  l1 l2 m1 m3 g Rho1 Rho2 A_1 A_2 E1 E2 I1 I2


% l1=0.1524;
% l2=0.3048;
% 
% % m2=1;
% m3=0;
% g=0;
% 
% A_1=0.00635^2/4*pi;
% A_2=A_1;
% Rho1=0.2534;
% Rho2=Rho1;
% E1=200e9;
% E2=E1;
% I1=7.9811E-11;
% I2=I1;
% m1=A_1*l1*Rho1;


% % l1=6;
% % l2=12;
% % g=0;
% % A_1=0.25^2/4*pi;
% % A_2=A_1;
% % Rho1=0.0007331;
% % Rho2=Rho1;      %similar %$
% % E1=30*10^6;
% % E2=E1;
% % I1=0.0001917476;
% % I2=I1;
% % m1=A_1*l1*Rho1;
% % m3=0;


% % 
% l1=0.15;
% l2=0.3;
% g=0;
% A_1=0.006^2/4*pi;
% A_2=A_1;
% Rho1=7.87e3;
% Rho2=Rho1;      %paper 748 E. M. BAKR and A. A. SHABANA 
% E1=0.2e12;
% E2=E1;
% I1=1e-10;
% I2=I1;
% m1=A_1*l1*Rho1;
% m3=m1;
% % 

l1=0.1524;
l2=0.304;

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


dt=0.01;
timeSpan=0:dt:1.6;
tsize=length(timeSpan);


% q_0=[pi/4;pi/8;0;0;0];
q_0=[0;0;0;0;0];
dq_0=[0;0;0;0;0];
mu_0=[0];

options=odeset('maxstep',1e-3);
%%
t0=clock;
[~,zIMM]=ode45(@IMM_dynamics,timeSpan,[q_0;dq_0;mu_0],options);
t1=clock;
timeIMM=etime(t1,t0);
disp(['IMM sim time: ' num2str(timeIMM) '(s)'])

% %%
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

IM.dq1=zIMM(:,6);
IM.dq2=zIMM(:,7);
IM.dq3=zIMM(:,8);
IM.dq4=zIMM(:,9);
IM.dq5=zIMM(:,10);


%
AM.q1=zAM(:,1);
AM.q2=zAM(:,2);
AM.q3=zAM(:,3);
AM.q4=zAM(:,4);
AM.q5=zAM(:,5);

AM.dq1=zAM(:,6);
AM.dq2=zAM(:,7);
AM.dq3=zAM(:,8);
AM.dq4=zAM(:,9);
AM.dq5=zAM(:,10);

%
EM.q1=zEM(:,1);
EM.q2=zEM(:,2);
EM.q3=zEM(:,3);
EM.q4=zEM(:,4);
EM.q5=zEM(:,5);


EM.dq1=zEM(:,6);
EM.dq2=zEM(:,7);
EM.dq3=zEM(:,8);
EM.dq4=zEM(:,9);
EM.dq5=zEM(:,10);

%
GM.q1=zGM(:,1);
GM.q2=zGM(:,2);
GM.q3=zGM(:,3);
GM.q4=zGM(:,4);
GM.q5=zGM(:,5);


GM.dq1=zGM(:,6);
GM.dq2=zGM(:,7);
GM.dq3=zGM(:,8);
GM.dq4=zGM(:,9);
GM.dq5=zGM(:,10);
%%
% 
Slider=cos(IM.q2).*(l2+IM.q4) + l1*cos(IM.q1);

% 
% % 
figure
plot(timeSpan,Slider,'color','#A2142F','linewidth',1.5)
grid minor
ylim([0.1 0.5])
%% 
%constraint
AE_IM = AEfunc(IM.dq1,IM.dq2,IM.dq4,l1,l2,IM.q1,IM.q2,IM.q4);
AE_AM = AEfunc(AM.dq1,AM.dq2,AM.dq4,l1,l2,AM.q1,AM.q2,AM.q4);
AE_EM = AEfunc(EM.dq1,EM.dq2,EM.dq4,l1,l2,EM.q1,EM.q2,EM.q4);
AE_GM = AEfunc(GM.dq1,GM.dq2,GM.dq4,l1,l2,GM.q1,GM.q2,GM.q4);

%%

figure
hold on; grid on
plot(timeSpan,AE_IM,'r-','linewidth',1.5)
plot(timeSpan,AE_AM,'g-','linewidth',1.5)
plot(timeSpan,AE_EM,'b-','linewidth',1.5)
plot(timeSpan,AE_GM,'k--','linewidth',1.5)
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


%%
% figure
% hold on; grid on
% plot(timeSpan,IM.q2,'r-','linewidth',2.5)
% plot(timeSpan,AM.q2,'g--','linewidth',2)
% plot(timeSpan,EM.q2,'b-','linewidth',1.5)
% plot(timeSpan,GM.q2,'k-','linewidth',1)
% legend('IM','AM','EM','GM')
% xlabel('time(s)')
% title('\delta_2')
% %
% %
% figure
% hold on; grid on
% plot(timeSpan,sin(IM.q5)*(l2/2),'r-','linewidth',2.5)
% plot(timeSpan,sin(AM.q5)*(l2/2),'g--','linewidth',2)
% plot(timeSpan,sin(EM.q5)*(l2/2),'b-','linewidth',1.5)
% plot(timeSpan,sin(GM.q5)*(l2/2),'k-','linewidth',1)
% legend('IM','AM','EM','GM')
% xlabel('time(s)')
% title('\delta_2')
%
%%
figure 
for i=1:length(timeSpan)
    
    XO(i)=0;
    YO(i)=0;
    
    XA(i)=XO(i)+l1*cos(IM.q1(i));
    YA(i)=YO(i)+l1*sin(IM.q1(i));
    
    XB(i)=XA(i)+(l2+IM.q4(i))*cos(IM.q2(i));
    YB(i)=YA(i)+(l2+IM.q4(i))*sin(IM.q2(i));
    
    
    P1x=[XO(i);XA(i)];
    P1y=[YO(i);YA(i)];
    
    P2x=[XA(i);XB(i)];
    P2y=[YA(i);YB(i)];
    
    
%     plot(xR,yR,'y-','linewidth',4)
%     hold on
%     plot(XO(i)+xd,YO(i)+yd,'r-','linewidth',4)
    
    plot(P1x,P1y,'b-','linewidth',1)
    hold on
    plot(P2x,P2y,'g-','linewidth',2)
    
   axis([-.5 1 -.5 .5])
    grid minor

    axis equal
    title(timeSpan(i))
    hold off
    pause(0.01)
    
end
    
