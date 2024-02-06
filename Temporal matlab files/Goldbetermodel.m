%%modelo de Goldbeter ritmo circadiano
vs=0.76;vm=0.65;vd=0.95;ks=0.38;k1=1.9;k2=1.3;
V1=3.2;V2=1.58;V3=5;V4=2.5;K1=1;K2=1;K3=1;K4=1;KI=1;
Km1=0.5;Kd=0.2;n=4;
%x=[x(1)=m,x(2)=p0,x(3)=P1,x(4)=p2,x(5)=pn]

F = @(t,x) [vs/(1+(x(5)/KI)^n) - vm*x(1)/(Km1+x(1));ks*x(1)-V1*x(2)/(K1+x(2))+V2*x(3)/(K2+x(3))...
    ;V1*x(2)/(K1+x(2))-V2*x(3)/(K2+x(3))-V3*x(3)/(K3+x(3))+V4*x(4)/(K4+x(4));V3*x(3)/(K3+x(3))-V4*x(4)/(K4+x(4))-k1*x(4)+k2*x(5)-vd*x(4)/(Kd+x(4))...
    ;k1*x(4)-k2*x(5)];
[t1,xa] = ode45(F,[0 200],[10 0 0 0 0]); % tiempo --> 0 a 1 , S0=5,E0=1,C0=0,P0=0

 figure(1)
plot(t1,xa(:,1));
hold on;
plot(t1,xa(:,2)+xa(:,3)+xa(:,4)+xa(:,5),'g');
hold on
% plot(t1,xa(:,3),'m');
% hold on
% plot(t1,xa(:,4),'y');
% hold on;
plot(t1,xa(:,5),'r');
% legend('[mRNA]','[PERprotein]','[PER-P]','[PERactiva]','PERnucleo','Location','best');
legend('[mRNA]','[PERtotal','PERnucleo','Location','best');
title('Goldbeter’s circadian oscillator model.');
xlabel('t(horas)'), ylabel('Concentracion');
grid on
axis([0 192 0 inf])
figure(2)
plot3(xa(:,1),xa(:,2)+xa(:,3)+xa(:,4)+xa(:,5),xa(:,5));
title('phase plane');
xlabel('mRNA'), ylabel('PERphos'),zlabel('PERnuclear');
grid on
