%Didier Gonze article based in goodwin model using Gillespie

NC=200;
     itime=100;
     timelimit=10000;
     tim=zeros(itime,1);
     mrn=zeros(itime,1);
     pro=zeros(itime,1);
     m1=zeros(itime,1);
     m2=zeros(itime,1);
     m3=zeros(itime,1);
     p1=zeros(itime,1);
     p2=zeros(itime,1);
     p3=zeros(itime,1);
%  parameters:reaction rates
v1=0.4; ka1=1.3; n=10; v2=0.6; k2=1; k3=0.7; v4=0.35; k4=1;
k5=0.7; v6=0.35; k6=1;
      
k=1; %counter for waitbar update 
alpha=10^2; %parameter for updating waitbar (increase for shorter runtime)
tic
%w=waitbar(0,'running SSA...');
for  is=1:NC
%  initialization
      M=5;
      P=0;
      I=0;
      %start time 
      time=0;
      izt=1;

     while time < timelimit
       %calculate rxn propensities
    a = [v1*(ka1^n)/(ka1^n+I^n) v2*M/(k2+M) k3*M v4*P/(k4+P) k5*P v6*I/(k6+I)];
    %combined rxn hazard
    a0 = sum(a);
    %calculate time to next event
    r1=rand;
    while r1 == 0,
        r1=rand;
    end
    tau=-(1./a0)*log(r1);
    %update time
    time = time + tau;
    %determine next reaction
    i=1; mu=0; amu=0; r2=rand;
    while amu < r2*a0,
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
    ztime=timelimit./itime;
       if time > izt*ztime
           time;
           izt;
           mrn(izt)=mrn(izt)+M;
           pro(izt)=pro(izt)+P;  
           if (is ==2)
              p1(izt)=P;
              m1(izt)=M;
                end
             if (is ==5)
              p2(izt)=P;
              m2(izt)=M;
                end
             if (is ==8)
              p3(izt)=P;
              m3(izt)=M;
             end
             izt=izt+1;
        end  
 %update waitbar
%     if t >= k*alpha*t_sample
%         k=k+1;
%         waitbar(t/t_end)
%     end     
     end
  %vectors containing steady state values for all cells
    hp(is)=P;
    hm(is)=M;
end
for jt=1:itime-1
       tim(jt)=ztime*jt;
        end
       prot=pro./NC;
       mrna=mrn./NC;
%close(w)
toc
figure(1)
        plot(tim,mrna)
        title('mRNA number average and 3 single cells')
        %legend('[P]','[M]','Location','best');
        xlabel('time')
        ylabel('concentration [nM]')
%         hold on 
%         plot(tim,m1,'r')
%         hold on
%         plot(tim,m2,'g')
%         hold on
%         plot(tim,m3,'m')
figure(2)
         plot(tim,prot)
        title('protein number: average and 3 single cells')
        xlabel('time')
        ylabel('concentration [nM]')
%         hold on 
%         plot(tim,p1,'r')
%         hold on
%         plot(tim,p2,'g')
%         hold on
%         plot(tim,p3,'m')
% figure(3)
%          hist(hp,10)
%         title('histogram of steady state protein number')
% 
% figure(4)
%         hist(hm,50)
%         title('histogram of steady state mRNA number')        