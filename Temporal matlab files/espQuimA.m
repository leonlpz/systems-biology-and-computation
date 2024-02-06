%modelo de produccion y degradacion de una especie quimica
 A0=1;
 k0=1;
% k1=3;
% a=(A0-k0/k1)*exp(-k1*t) + k0/k1;
% b=log(A0*k0/k1 - (k0/k1)^2) - k1*t;
 c=1./(-2*k0*t+1/A0);
 t=0:0.1:10;
% figure;
% plot(t,a);
% figure;
 plot(t,c);
%-------------------------------------------------
%modelo de interaccion entre las especies quimicas A y B
% k3=12;
% k2=2;
% k1=9;
% %x=[x(1)=A,x(2)=B]
% f = @(t,x) [k3*x(2)-k1*x(1);k1*x(1)-(k3+k2)*x(2)];
% [t,xa] = ode45(f,[0 3],[0 10]); % tiempo --> 0:t , A0=0, B0=10
% t1=1/(k1+k3)
% t2=1/k2
% %subplot(2,1,1);
% figure;
% plot(t,xa(:,1));
% hold on;
% plot(t,xa(:,2),'color','r');
% legend('[A]','[B]','Location','best');
% title('Concentraciones A y B');
% xlabel('t'), ylabel('[A],[B]');
% grid on
% %subplot(2,1,2);
% figure;
% plot(xa(:,1),xa(:,2));
% title('plano fase A y B');
% xlabel('A'), ylabel('B');