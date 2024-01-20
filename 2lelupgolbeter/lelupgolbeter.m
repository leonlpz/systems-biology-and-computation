%Lelup-Golbeter --> 16 variable circadian rhythm model
%mRNA per,cry,Bmal1---------------------------------------------
%Mp=X(1), Mc=x(2), Mb=x(3)
%per, cry phosphorilated no-phosphorilated----------------------
%Pc=x(4), Cc=x(5), Pcpx(6), Ccp=x(7)
%complex PER-CRY phosphorilated no-phosphorilated---------------
%PCc=x(8), PCn=x(9), PCcp=x(10), PCnp=x(11)
%Bmal1 phosphorilated no-phosphorilated citosol nucleus---------
%Bc=x(12), Bcp=x(13), Bn=x(14), Bnp=x(15)
%PER-CRY CLOCK-BMAL1 inactive complex nucleus-------------------
%In=x(16)
k1=0.4;k2=0.2;k3=0.4;k4=0.2;k5=0.4;k6=0.2;k7=0.5;k8=0.1;
KAP=0.7;KAC=0.6;KIB=2;kdmc=0.01;kdmb=0.01;kdmp=0.01;kdn=0.1;kdnc=0.12;kd=0.3;
kdp=0.1;kp=0.1;kmb=0.9;kmc=0.4;kmp=0.31;ksb=0.12;ksc=1.6;ksp=0.6;
m=2;n=4;
v1b=0.5;v1c=0.6;v1p=0.4;v1pc=0.4;v2b=0.1;v2c=0.1;v2p=0.3;v2pc=0.1;
v3b=0.5;v3pc=0.4;v4b=0.2;v4pc=0.1;vp=0.4;vdbc=0.5;vdbn=0.6;vdcc=0.7;
vdin=0.8;vdpc=0.7;vdpcc=0.7;vdpcn=0.7;
vmb=0.8;vmc=1.0;vmp=1.1;vsb=1.0;vsc=1.1;vsp=1.5;

F = @(t,x) [vsp*(x(14))^n/((KAP)^n+x(14)^n)-vmp*x(1)/(kmp+x(1))-kdmp*x(1)...[1]
           ;vsc*(x(14))^n/((KAC)^n+x(14)^n)-vmc*x(2)/(kmc+x(2))-kdmc*x(2)...[2]
           ;vsb*(KIB^m)/((KIB^m)+x(14)^m)-vmb*x(3)/(kmb + x(3))-kdmb*x(3)...[3]
           ;ksp*x(1)-v1p*x(4)/(kp+x(4))+v2p*x(6)/(kdp+x(6))-k3*x(4)*x(5)+k4*x(8)-kdn*x(4)...[4]
           ;ksc*x(2)-v1c*x(5)/(kp+x(5))-v2c*x(7)/(kdp+x(7))-k3*x(4)*x(5)+k4*x(8)-kdnc*x(5)...[5]
           ;v1p*x(4)/(kp+x(4))-v2p*x(6)/(kdp+x(6))-vdpc*x(6)/(kd+x(6))-kdn*x(6)...[6]
           ;v1c*x(5)/(kp+x(5))-v2c*x(7)/(kdp+x(7))-vdcc*x(7)/(kd+x(7))-kdn*x(7)...[7]
           ;-v1pc*x(8)/(kp+x(8))+v2c*x(10)/(kdp+x(10))+k3*x(4)*x(5)-(k4+k1+kdn)*x(8)+k2*x(9)...[8]
           ;-v3pc*x(9)/(kp+x(9))+v4pc*x(11)/(kdp+x(11))+k1*x(8)-(k2+kdn)*x(9)-k7*x(14)*x(9)+k8*x(16)...[9]
           ;v1pc*x(8)/(kp+x(8))-v2pc*x(10)/(kdp+x(10))-vdpcc*x(10)/(kd+x(10))-kdn*x(10)...[10]
           ;v3pc*x(9)/(kp+x(9))-v4pc*x(11)/(kdp+x(11))-vdpcn*x(11)/(kd+x(11))-kdn*x(11)...[11]
           ;ksb*x(3)-v1b*x(12)/(kp+x(12))+v2b*x(13)/(kdp+x(13))-(k5+kdn)*x(12)+k6*x(14)...[12]
           ;v1b*x(12)/(kp+x(12))-v2b*x(13)/(kdp+x(13))-vdbc*x(13)/(kd+x(13))-kdn*x(13)...[13]           
           ;-v3b*x(14)/(kp+x(14))+v4b*x(15)/(kdp+x(15))+k5*x(12)-(k6+kdn)*x(14)-k7*x(14)*x(9)+k8*x(16)...[14]
           ;v3b*x(13)/(kp+x(13))-v4b*x(15)/(kdp+x(15))-vdbn*x(15)/(kd+x(15))-kdn*x(15)...[15]
           ;-k8*x(16)+k7*x(14)*x(9)-vdin*x(16)/(kd+x(16))-kdn*x(16)];%[16]
[t,xa] = ode45(F,[0 500],[10 10 10 0 0 0 0 0 0 0 0 0 0 0 0 0]);

figure(1);
plot(t,xa(:,1),'r');
hold on;
plot(t,xa(:,4)+xa(:,6)+xa(:,8)+xa(:,9)+xa(:,10)+xa(:,11));
hold on;
plot(t,xa(:,8)+xa(:,9)+xa(:,10)+xa(:,11),'g');
hold on;
legend('permRNA','totalPER','PERnucleo','Location','best');
%title('Lelup-Golbeter mo');
xlabel('t(horas)'), ylabel('Concentraciones [nM]');
grid on
axis([0 192 0 inf]);
figure(2);
plot3(xa(:,1),xa(:,4)+xa(:,6)+xa(:,8)+xa(:,9)+xa(:,10)+xa(:,11),+xa(:,6)+xa(:,8)+xa(:,9)+xa(:,10)+xa(:,11));
hold on;
title('pplane mRNA PER ');
xlabel('permRNA'), ylabel('totalPER'),zlabel('nuclearPER');
grid on
figure(3)
plot(t,xa(:,1),'r');
hold on;
plot(t,xa(:,2));
hold on;
plot(t,xa(:,3),'g');
hold on;
legend('permRNA','crymRNA','BmalmRNA','Location','best');
%title('Lelup-Golbeter mo');
xlabel('t(horas)'), ylabel('Concentraciones [nM]');
grid on
axis([0 192 0 inf]);
figure(4)
plot3(xa(:,1),xa(:,2),xa(:,3));
hold on;
title('pplane per,cry,Bmal mRNA´s ');
xlabel('per'), ylabel('cry'),zlabel('Bmal');
grid on
figure(5)
plot(t,xa(:,4)+xa(:,6)+xa(:,8)+xa(:,9)+xa(:,10)+xa(:,11),'r');
hold on;
plot(t,xa(:,5)+xa(:,7)+xa(:,8)+xa(:,9)+xa(:,10)+xa(:,11));
hold on;
plot(t,xa(:,12)+xa(:,13)+xa(:,14)+xa(:,15),'g');
hold on;
legend('PER','CRY','BMAL','Location','best');
%title('Lelup-Golbeter mo');
xlabel('t(horas)'), ylabel('Concentraciones [nM]');
grid on
axis([0 192 0 inf]);
figure(6)
plot3(xa(:,4)+xa(:,6)+xa(:,8)+xa(:,9)+xa(:,10)+xa(:,11),xa(:,5)+xa(:,7)+xa(:,8)+xa(:,9)+xa(:,10)+xa(:,11)...
    ,xa(:,12)+xa(:,13)+xa(:,14)+xa(:,15));
hold on;
title('pplane PER,CRY,BMAL protein ');
xlabel('PER'), ylabel('CRY'),zlabel('BMAL');
grid on
