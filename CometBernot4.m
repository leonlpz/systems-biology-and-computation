%Comet-Bernot --> reduced 4 variable model circadian rhythm
%protein per-cito Pc=x(1); protein cry-cito Cc=x(2); complex cito PCc=x(3);
%complex nucl PCn=x(4)
K=0.4; n=30; v1=2; v2=2.2; k1=0.08; k2=0.06;
k3=0.08; k4=0.06; kd1=0.05; kd2=0.05;kd3=0.05;
kd4=0.25;
F = @(t,x) [v1*(K^n)/((K^n)+x(4)^n)-k3*x(1)*x(2)+k4*x(3)-kd1*x(1);v2*(K^n)/((K^n)+x(4)^n)-k3*x(1)*x(2)+k4*x(3)-kd2*x(2)...
    ;k3*x(1)*x(2)-(k4+k1+kd3)*x(3)+k2*x(4);k1*x(3)-(k2+kd4)*x(4)];
[t,xa] = ode45(F,[0 192],[0 0 0 0]);

figure(1);
plot(t,xa(:,1));
hold on;
plot(t,xa(:,2),'r');
hold on;
plot(t,xa(:,3),'g');
hold on;
plot(t,xa(:,4),'m');
hold on;
legend('PERc','CRYc','PER-CRY_c','PER-CRY_n','Location','best');
title('PER-CRY proteins');
xlabel('time(h)'), ylabel('Concentraciones [nM]');
grid on
figure(2);
plot(xa(:,1),xa(:,2));
hold on;
title('pplane proteins');
xlabel('PERc'), ylabel('CRYc');
grid on
figure(3)
plot3(xa(:,1),xa(:,2),xa(:,3));
hold on;
title('pplane proteins');
xlabel('PERc'), ylabel('CRYc'),zlabel('PER-CRYc');
grid on
figure(4)
plot3(xa(:,1),xa(:,2),xa(:,4));
hold on;
title('pplane proteins');
xlabel('PERc'), ylabel('CRYc'),zlabel('PER-CRYn');
grid on
