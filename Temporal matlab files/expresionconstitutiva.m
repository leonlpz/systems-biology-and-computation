%modelo de expresión constitutiva 
%mRNA'=k1[ADN]-k2mRNA
%protein'=k3*mRNA-k4*protein
%realimentacion negativa tipo Michaelis Menten:
%mRNA'=k1*(k/k+protein)[ADN]-k2mRNA
%protein'=k3*m-k4*protein
%[ADN]=1
k=1000;%controla la afinidad del lazo realimentado 
k1=10000;
k2=10;
k3=3;
k4=2;
F = @(t,x) [k1-k2*x(1);k3*x(1)-k4*x(2)];
[t1,xa] = ode45(F,[0 2],[1 1]);

G = @(t,x) [k1*(k/(k+x(2)))-k2*x(1);k3*x(1)-k4*x(2)];
[t2,xb] = ode45(G,[0 2],[1 1]);

 figure(1);
 plot(t1,xa(:,1));
 hold on;
 plot(t1,xa(:,2),'r');
 hold on;
 plot(t2,xb(:,1),'--b');
 hold on;
 plot(t2,xb(:,2),'--r');
 legend('[mRNA]','[proteina]','[mRNA reali(-)]','[proteina reali(-)]','Location','best');
 title('Expresion constitutiva');
 xlabel('tiempo'), ylabel('Concentracion');
 grid on
 
 figure(2);
 plot(xa(:,1),xa(:,2),'b');
 hold on;
 plot(xb(:,1),xb(:,2),'r');
 hold on;
 clear all
 clc