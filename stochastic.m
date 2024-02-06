%Didier Gonze article based in goodwin model using Gillespie

NC=200;%cell number
     t_end=1000; %end time
     t_sample=1; %sample interval for gathering data
%  parameters:reaction rates
% ka1,k2,k4 and k6 are average disociation constants
v1=1; ka1=1.3; n=4; v2=0.6; k2=1; k3=0.7; v4=1; k4=1;
k5=0.7; v6=0.35; k6=1; OM=1000;%volumen celula
tic
for  is=1:1:NC
%  initialization
      M=0;
      P=0;
      I=0;
      t=0; %start time 
      %%%arrays to store results
j=1; %counter for arrays
t_array(1,t_end/t_sample)=zeros(1); t_array(1,j)=t; %time array and initial value
M_array(1,t_end/t_sample)=zeros(1); M_array(1,j)=M; %substrate array and initial value
P_array(1,t_end/t_sample)=zeros(1); P_array(1,j)=P; %product array and initial value
I_array(1,t_end/t_sample)=zeros(1); I_array(1,j)=I;
%      mrna=zeros(length(M_array),1);
%      pro=zeros(length(P_array),1);
%      time=zeros(length(t_array),1);

     while t < t_end
       %calculate rxn propensities
    a = [(v1*OM)*((ka1*OM)^n)/((ka1*OM)^n+I^n) v2*OM*M/(k2+M) k3*M v4*OM*P/(k4*OM+P) k5*P v6*OM*I/(k6*OM+I)];
    %combined rxn hazard
    a0 = sum(a);
    %calculate time to next event
    r1=rand;
    while r1 == 0
        r1=rand;
    end
    tau=-(1./a0)*log(r1);
    %update time
    t = t + tau;
    %determine next reaction
    i=1; mu=0; amu=0; r2=rand;
    while amu < r2*a0
        mu = mu + 1;
        amu = amu + a(i); 
        i = i + 1;				
    end   
      if (mu == 1) 
      M=M+1;
      elseif(mu == 2)
      M=M-1;
      elseif(mu == 3) 
      P=P+1;
      elseif(mu == 4) 
      P=P-1;
      elseif(mu == 5) 
      I=I+1;
      elseif(mu == 6) 
      I=I-1;
      end 
      %store/output time and species
    if t >= j*t_sample
        j=j+1;
        t_array(1,j)=j;
        M_array(1,j)=M;
        P_array(1,j)=P;
        I_array(1,j)=I;
    end 
    
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
aa=round(NC*rand(1));
bb=round(NC*rand(1));
figure() 
plot(time(:,aa),pro(:,aa))
hold on
plot(time(:,aa),mrna(:,aa),'r')
        title('mRNA y proteina 1 rand')
        legend('[P]','[M]','Location','best');
        xlabel('time')
        ylabel('Numero de moleculas')
        axis([0 192 0 inf])

figure()
plot(mrna(:,aa),pro(:,aa))
title('pplane 1 rand mRNA - proteina');
xlabel('mRNA'), ylabel('protein');
grid on 
figure()
plot3(mrna(:,aa),pro(:,aa),inh(:,aa))
title('pplane mRNA - proteins');
xlabel('mRNA'), ylabel('protein'),zlabel('inhibitor complex');
grid on
figure() 
plot(time(:,bb),pro(:,bb))
hold on
plot(time(:,bb),mrna(:,bb),'r')
        title('mRNA y proteina 1 rand')
        legend('[P]','[M]','Location','best');
        xlabel('time')
        ylabel('Numero de moleculas')
        axis([0 192 0 inf])
        
figure()
plot(mrna(:,bb),pro(:,bb))
title('pplane 1 rand mRNA - protein');
xlabel('mRNA'), ylabel('protein');
grid on
%axis([0 14 0 30])        
figure()
plot3(mrna(:,bb),pro(:,bb),inh(:,bb))
title('pplane mRNA - proteins');
xlabel('mRNA'), ylabel('protein'),zlabel('inhibitor complex');
grid on
%figure()
% plot(tiempo,Iprom,'g')
%         title('complex')
%         %legend('[P]','[M]','Location','best');
%         xlabel('time')
%         ylabel('Number of molecules')
%         axis([0 192 0 inf])        
figure()
plot(tiempo,Pprom)
hold on
plot(tiempo,Mprom,'r')
        title('mRNA proteina promedio')
        legend('[P]','[M]','Location','best');
        xlabel('time')
        ylabel('Numero de moleculas')
        axis([0 192 0 inf])
figure()
plot(Mprom,Pprom)
title('pplane mRNA - protein promedio');
xlabel('mRNA'), ylabel('Protein');
grid on
figure()
plot3(Mprom,Pprom,Iprom)
title('pplane mRNA - protein -complex promedio');
xlabel('mRNA'), ylabel('Protein'),zlabel('Inhibition complex');
grid on

%---------------------------------------------------------        
figure()
[picos ,locs]=findpeaks(mrna(:,aa),'minpeakdistance',15);
peakInterval = diff(time(locs,aa));
hist(peakInterval)
grid on
xlabel('time Intervals')
ylabel('Occurrence')
title('Histogram of Peak Intervals (h)')
%---------------------------------------------------------
%-------------------------------------------------------        
% figure()
% [pks, locs2] = findpeaks(mrna(:,bb),time(:,bb),'MinPeakProminence',40);
% peakInterval = diff(locs2);
% hist(peakInterval)
% grid on
% xlabel('time Intervals')
% ylabel('Occurrence')
% title('Histogram of Peak Intervals (h)')
%---------------------------------------------------------
%AverageDistance_Peaks1 = mean(diff(locs));
% AverageDistance_Peaks2 = mean(diff(locs2));