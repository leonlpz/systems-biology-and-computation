%stochastic model for Drosophila circadian rhythms
NC=200;%cell number
     t_end=1000; %end time
     t_sample=1; %sample interval for gathering data
%  parameters:reaction rates
 n=4;vs=1;vm=0.7;vdp=2;vdt=2;ks=0.9;k1=0.6;k2=0.2;k3=1.2;k4=0.6;
 km=0.2;ki=1;kdpt=0.2;knpt=2;v1=8;v2=1;v3=8;v4=1;kdcn=0.01;
OM=1;%volumen celula
tic
for  is=1:NC
%  initialization
      Mp=10;
      P0=0;
      P1=0;
      P2=0;
      Mt=10;
      T0=0;
      T1=0;
      T2=0;
      C=0;
      Cn=0;
      t=0; %start time
      %%%arrays to store results
j=1; %counter for arrays
t_array(1,t_end/t_sample)=zeros(1); t_array(1,j)=t; %time array and initial value
M_array(1,t_end/t_sample)=zeros(1); M_array(1,j)=Mp; %substrate array and initial value
P_array(1,t_end/t_sample)=zeros(1); P_array(1,j)=Cn; %product array and initial value
I_array(1,t_end/t_sample)=zeros(1); I_array(1,j)=Mt;
%      mrna=zeros(length(M_array),1);
%      pro=zeros(length(P_array),1);
%      time=zeros(length(t_array),1);

     while t < t_end
       %calculate rxn propensities
w1=vs*OM*(ki*OM)^n/((ki*OM)^n+Cn^n);w2=vm*OM*Mp/(km*OM+Mp);w3=ks*Mp;w4=v1*OM*P0/(knpt*OM+P0);w5=v2*OM*P1/(knpt*OM+P1);
w6=v3*OM*P1/(knpt*OM+P1);w7=v4*OM*P2/(knpt*OM+P2);w8=k3*P2*T2/OM;w9=k4*C;w10=vdp*OM*P2/(kdpt*OM+P2);
w11=vs*OM*(ki*OM)^n/((ki*OM)^n+Cn^n);w12=vm*OM*Mt/(km*OM+Mt);w13=ks*Mt;w14=v1*OM*T0/(knpt*OM+T0);w15=v2*OM*T1/(knpt*OM+T1);
w16=v3*OM*T1/(knpt*OM+T1);w17=v4*OM*T2/(knpt*OM+T2);w18=vdt*OM*T2/(kdpt*OM+T2);w19=k1*C;w20=k2*Cn;
w21=kdcn*Mp;w22=kdcn*P0;w23=kdcn*P1;w24=kdcn*P2;w25=kdcn*Mt;
w26=kdcn*T0;w27=kdcn*T1;w28=kdcn*T2;w29=kdcn*C;w30=kdcn*Cn;
    a = [w1 w2 w3 w4 w5 w6 w7 w8 w9 w10 w11 w12 w13 w14 w15 w16 w17 w18 w19 w20 w21 w22 w23 w24 w25 w26 w27 w28 w29 w30];
    %combined rxn hazard
    a0 = sum(a);
    %calculate time to next event
    r1=rand;
    while r1 == 0,
        r1=rand;
    end
    tau=-(1./a0)*log(r1);
    %update time
    t = t + tau;
    %determine next reaction
    i=1; mu=0; amu=0; r2=rand;
    while amu < r2*a0,
        mu = mu + 1;
        amu = amu + a(i); 
        i = i + 1;				
    end   
      if (mu == 1) 
      Mp=Mp+1;
      elseif(mu == 2)
      Mp=Mp-1;
      elseif(mu == 3) 
      Mp=Mp-1;
      P0=P0+1;
      elseif(mu == 4) 
      P0=P0-1;
      P1=P1+1;
      elseif(mu == 5) 
      P1=P1-1;
      P0=P0+1;
      elseif(mu == 6) 
      P1=P1-1;
      P2=P2+1;
      elseif(mu == 7)
      P2=P2-1;
      P1=P1+1;
      elseif(mu == 8) 
      P2=P2-1;
      T2=T2-1;
      C=C+1;
      elseif(mu == 9) 
      C=C-1;
      P2=P2+1;
      T2=T2+1;
      elseif(mu == 10) 
      P2=P2-1;
      elseif(mu == 11) 
      Mt=Mt+1;
      elseif(mu == 12)
      Mt=Mt-1;
      elseif(mu == 13) 
      Mt=Mt-1;
      T0=T0+1;
      elseif(mu == 14) 
      T0=T0-1;
      T1=T1+1;
      elseif(mu == 15) 
      T1=T1-1;
      T0=T0+1;
      elseif(mu == 16) 
      T1=T1-1;
      T2=T2+1;
      elseif(mu == 17)
      T2=T2-1;
      T1=T1+1;
      elseif(mu == 18) 
      T2=T2-1;
      elseif(mu == 19) 
      C=C-1;
      Cn=Cn+1;
      elseif(mu == 20) 
      Cn=Cn-1;
      C=C+1;
      elseif(mu == 21) 
      Mp=Mp-1;
      elseif(mu == 22)
      P0=P0-1;
      elseif(mu == 23) 
      P1=P1-1;
      elseif(mu == 24) 
      P2=P2-1;
      elseif(mu == 25) 
      Mt=Mt-1;
      elseif(mu == 26) 
      T0=T0-1;
      elseif(mu == 27)
      T1=T1-1;
      elseif(mu == 28) 
      T2=T2-1;
      elseif(mu == 29) 
      C=C-1;
      elseif(mu == 30) 
      Cn=Cn-1;
      end 
      %store/output time and species
    if t >= j*t_sample
        j=j+1;
        t_array(1,j)=j;
        M_array(1,j)=Mp;
        P_array(1,j)=Cn;
        I_array(1,j)=Mt;
    end 
    %%%%
     end 
    mrna(:,is)=M_array;
    pro(:,is)=P_array;
    time(:,is)=t_array;
    inh(:,is)=I_array;
          
end
     Mprom=mean(mrna');
     Pprom=mean(pro');
     tiempo=mean(time');
     Iprom=mean(inh');
toc;
if toc > 60
    X=toc/60;
    elapsed = ['Elapsed time in minutes: ',num2str(X)];
    disp(elapsed)
end
% Nsamps = t_end;
% fsamp = 1000;
% Tsamp = 1/fsamp;
% tt = (0:Nsamps-1)*Tsamp;
% % Choose FFT size and calculate spectrum
% Nfft = 1000;%1024;
% [Pxx,f] = pwelch(Pprom,gausswin(Nfft),Nfft/2,Nfft,fsamp);

figure(1) 
plot(t_array,P_array)
 title('per-tim nuclear complex')
        %legend('[P]','[M]','Location','best');
        xlabel('time')
        ylabel('Number of molecules')
        axis([0 192 0 inf])
figure(2)
plot(t_array,M_array,'r')
        title('per mRNA')
        %legend('[P]','[M]','Location','best');
        xlabel('time')
        ylabel('Number of molecules')
        axis([0 192 0 inf])
figure(3)
plot(tiempo,Iprom,'g')
        title('tim mRNA average over cells')
        %legend('[P]','[M]','Location','best');
        xlabel('time')
        ylabel('Number of molecules')
        axis([0 192 0 inf])        
figure(4)
plot(tiempo,Pprom)
 title('per-tim nuclear complex averaged over cells')
        %legend('[P]','[M]','Location','best');
        xlabel('time')
        ylabel('Number of molecules')
        axis([0 192 0 inf])
figure(5)
plot(tiempo,Mprom,'r')
        title('per mRNA averaged over cells')
        %legend('[P]','[M]','Location','best');
        xlabel('time')
        ylabel('Number of molecules')
        axis([0 192 0 inf])
figure(6)
plot(Mprom,Pprom)
title('pplane mRNA - proteins');
xlabel('per mRNA'), ylabel('per-tim nuclear complex');
grid on
%axis([0 14 0 30])
figure(7)
plot3(Mprom,Iprom,Pprom)
title('pplane mRNA - proteins averaged over cells');
xlabel('per mRNA'), ylabel('tim mRNA'),zlabel('per-tim nuclear complex');
grid on
%axis([0 5 0 10])
% figure(7)
% hist(Mprom,50)
% figure(8)
% hist(Pprom,50)