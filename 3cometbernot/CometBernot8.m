%Comet-Bernot --> 8 variable reduced circadian rhythm model
%per_cito Pc=x(1); cry_cito Cc=x(2); complex_cito PCc=x(3);
%complex_nucl PCn=x(4);per_cito-phos Pcp=x(5); cry_cito-phos Ccp=x(6);
%complex_cito-phos;PCcp=x(7);complex_nuc-phos PCnp=x(8)
K=0.01; n=10; v1=1.7; v2=2.5; k1=0.08; k2=0.06;
k3=0.08; k4=0.06; kd1=0.05; kd2=0.05;kd3=0.05;
kd4=0.30;v1p=0.01;v1c=0.01;v2p=0.3;v2c=0.1;v1pc=0.4;
v2pc=0.1;v3pc=0.4;v4pc=0.2;vdpc=0.7;kdn=0.01;
vdcc=0.7;vdpcc=0.7;vdpcn=0.7;kd=0.3;kp=5;kdp=5;

F = @(t,x) [v1*(K^n)/((K^n)+x(4)^n)-v1p*x(1)/(kp+x(1))+v2p*x(5)/(kdp+x(5))-k3*x(1)*x(2)+k4*x(3)-kd1*x(1)...(1)
    ;v2*(K^n)/((K^n)+x(4)^n)-v1c*x(2)/(kp+x(2))+v2c*x(6)/(kdp+x(6))-k3*x(1)*x(2)+k4*x(3)-kd2*x(2)...(2)
    ;-v1pc*x(3)/(kp+x(4))+v2c*x(7)/(kdp+x(7))+k3*x(1)*x(2)-(k4+k1+kd3)*x(3)+k2*x(4)...(3)
    ;-v3pc*x(4)/(kp+x(4))+v4pc*x(8)/(kdp+x(8))+k1*x(3)-(k2+kd4)*x(4)...(4)
    ;v1p*x(1)/(kp+x(1))-v2p*x(5)/(kdp+x(5))-vdpc*x(5)/(kd+x(5))-kdn*x(5)...(5)
    ;v1c*x(2)/(kp+x(2))-v2c*x(6)/(kdp+x(6))-vdcc*x(6)/(kd+x(6))-kdn*x(6)...(6)
    ;v1pc*x(3)/(kp+x(3))-v2pc*x(7)/(kdp+x(7))-vdpcc*x(7)/(kd+x(7))-kdn*x(7)...(7)
    ;v3pc*x(4)/(kp+x(4))-v4pc*x(8)/(kdp+x(8))-vdpcn*x(8)/(kd+x(8))-kdn*x(8)];%(8)
[t,xa] = ode45(F,[0 192],[1 1 1 1 0 0 0 0]);

figure(1);
plot(t,xa(:,1));
hold on;
plot(t,xa(:,2),'r');
hold on;
legend('PERc','CRYc','Location','best');
title('PER CRY proteins');
xlabel('time(h)'), ylabel('Concentraciones [nM]');
grid on
figure(2);
plot(xa(:,1),xa(:,2));
hold on;
title('pplane PER CRY proteins');
xlabel('PERc'), ylabel('CRYc');
grid on
