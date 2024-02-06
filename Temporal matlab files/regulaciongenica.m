%%x'==X(1) , y'==x(2)
F = @(t,x) [((0.01+0.72*x(1)^3)/(1+8*x(1)^3+x(2)^3))-0.1*x(1);((0.01+0.72*x(2)^3)/(1+8*x(2)^3+x(1)^3))-0.1*x(2)];
[t,xa] = ode45(F,[0 500],[2 1]);

figure(1);
plot(t,xa(:,1));
hold on;
plot(t,xa(:,2),'r');
hold on
legend('[x=G1]','[y=G2]','Location','best');
%title('');
xlabel('time'), ylabel('Concentración de mRNA');
grid on

figure(2);
plot(xa(:,1),xa(:,2));
hold on;
%legend('[x=G1]','[y=G2]','Location','best');
title('Plano fase');
xlabel('Concentración mRNA G1'), ylabel('Concentración mRNA G2');
grid on



