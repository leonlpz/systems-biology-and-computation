%Didier Gonze ---> Goodwin model
%parameters:
%mRNA, protein, inhibitor, of "per" or "cry" molecule
%mRNA=x(1) ; Protein=x(2) ; Inhibitor=x(3)
v1=0.4; K1=1.3; n=10; v2=0.6; K2=1; k3=0.7; v4=0.35; K4=1;
k5=0.7; v6=0.35; K6=1;
time=500;
F = @(t,x) [v1*(K1^n)/((K1^n)+x(3)^n) - v2*x(1)/(K2+x(1))...
           ;k3*x(1)-v4*x(2)/(K4+x(2));k5*x(2)-v6*x(3)/(K6+x(3))];
[t,xa] = ode45(F,[0 time],[0 0 0]);

figure(1);
plot(t,xa(:,1),'r');
hold on;
plot(t,xa(:,2));
hold on;
% plot(t,xa(:,3),'g');
% hold on;
legend('mRNA','Protein','Location','best');
title('Concentracion [nM] "per"');
xlabel('time(h)'), ylabel('[M],[P]');
grid on
axis([0 196 0 inf]);
figure(2);
plot3(xa(:,1),xa(:,2),xa(:,3));
hold on;
title('pplane mRNA proteins');
xlabel('mRNA'), ylabel('Protein'),zlabel('inhibitor complex');
grid on
figure(3);
plot(t,xa(:,3));
xlabel('time(h)'), ylabel('inhibitor complex');
axis([0 196 0 inf]);
grid on
 figure(4)
 plot(xa(:,1),xa(:,2));
 xlabel('mRNA'), ylabel('Protein')
 hold on;
% hist(xa(:,1),50)
% title('mRNA');
% figure(5)
% hist(xa(:,2),50)
% title('PERc');
%FFT
% L=length(t);
% Fs=100;
% Y = fft(xa(:,1));
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% %P1(2:end-1) = 2*P1(2:end-1);
% f = Fs*(0:(L/2))/L;
% figure(6)
% plot(f,P1) 
% title('Single-Sided Amplitude Spectrum of mRNA')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')