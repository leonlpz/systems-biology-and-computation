%%modelo de cinetica quimica
k11=1;%k11=k-1
k2=10;
k1=30;
k22=20;%k22=k-2
et=1;
Vmax=k2*et;
km=(k11+k2)/k1;
vf=k2*et;
vr=k11*et;
num1=(k1/(k11+k2))*vf;
num2=(k22/(k11+k2))*vr;
den1=k1/(k1+k2);
den2=k22/(k11+k2);
%t=linespace(0:0.1)
%x=[x(1)=S,x(2)=C,x(3)=P]
%E=et-x(2)
F = @(t,x) [k11*x(2)-k1*x(1)*(et-x(2));k1*x(1)*(et-x(2))-(k11+k2)*x(2);k2*x(2)];
[t1,xa] = ode45(F,[0 2],[5 0 0]); % tiempo --> 0 a 1 , S0=5,E0=1,C0=0,P0=0

f=@(t,x) [-Vmax*x(1)/(km+x(1));Vmax*x(1)/(km+x(1))];
[t2,xb] = ode45(f,[0 2],[5 0]);

G=@(t,x) [k11*x(2)-k1*x(1)*(et-x(2));k1*x(1)*(et-x(2))+k22*x(3)*(et-x(2))-(k11+k2)*x(2);k2*x(2)-k22*x(3)*(et-x(2))];
[t3,xc]=ode45(G,[0 2],[5 0 0]);

g=@(t,x) [-(num1*x(1)-num2*x(2))/(1+den1*x(1)+den2*x(2));(num1*x(1)-num2*x(2))/(1+den1*x(1)+den2*x(2))];
[t4,xd]=ode45(g,[0 2],[5 0]);

% figure(1);
% plot(t1,xa(:,1));
% hold on;
% plot(t1,xa(:,2),'r');
% hold on
% plot(t1,xa(:,3),'g');
% hold on
% plot(t1,et-xa(:,2),'m');
% legend('[S]','[C]','[P]','[E]','Location','best');
% title('Concentraciones S, E, C y P');
% xlabel('time'), ylabel('[S],[C],[P],[E]');
% grid on

% figure(2);
% plot(t1,xa(:,1));
% hold on;
% plot(t2,xb(:,1),'--b');
% hold on;
% plot(t1,xa(:,3),'r');
% hold on;
% plot(t2,xb(:,2),':r');
% legend('[S]full','[S]reduced','[P]full','[P]reduced','Location','best');
% title('Concentraciones S-->P');
% xlabel('time'), ylabel('[S],[P]');
% grid on

% figure(3);
% plot(t1,xa(:,1));
% hold on;
% plot(t3,xc(:,1),'--b');
% hold on;
% plot(t1,xa(:,3),'r');
% hold on;
% plot(t3,xc(:,3),'--r');
% hold on;
% legend('[Snr]','[Sr]','[Pnr]','[Pr]','Location','best');
% title('Concentraciones Snorev, Srev, Pnorev y Prev');
% xlabel('time'), ylabel('[Snr],[Sr],[Pnr],[Pr]');
% grid on

figure(4);
plot(t1,xa(:,1));
hold on;
plot(t2,xb(:,1),'--b');
hold on;
plot(t3,xc(:,1),'*b');
hold on;
plot(t4,xd(:,1),'ob');
hold on;
plot(t1,xa(:,3),'r');
hold on;
plot(t2,xb(:,2),'--r');
hold on;
plot(t3,xc(:,3),'*r');
hold on;
plot(t4,xd(:,2),'or');
hold on;
legend('[S]full','[S]reduced','[S]fullrever','[S]redurever','[P]full','[P]reduced','[P]fullrever','[P]redurever','Location','best');
title('Concentraciones S-->P');
xlabel('time'), ylabel('[S],[P]');
grid on
