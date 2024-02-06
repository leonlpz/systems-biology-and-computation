%stochastic implementation (Gillespie SSA) of enzimatic reactions
%Figure 3.3 Ingalls
%the network#1 is simply  S+E <--> C --> E+P . with propensity: a1=k1*N_s*N_e,
%a_1=k_1*N_c, b=k2*N_c
%where N_x is the
%number of molecules of X
%Network#2  S --> P , propensity: a=k*N_s
% reaction s+e=c (c1) % k1=30;
% reaction c=s+e (c2) % k11=1;%k11=k-1
% reaction c=e+p (c3) % k2=10;
% parameters of the reaction system, nens=nb of cells, timelimit= length
%    of running time, itime is length of file which gives time 
% dependence of proteins or averages
clc; clear all; close all;
 
     NC=1;
     t_end=1000; %end time
     t_sample=1; %sample interval for gathering data
%  parameters:reaction rates
      k1=30;
      k2=1;
      k3=10;
k=1; %counter for waitbar update 
alpha=10^2; %parameter for updating waitbar (increase for shorter runtime)
tic
w=waitbar(0,'running SSA...');
for  is=1:NC
%  initialization
      s=5;
      e=1;
      c=0;
      p=0;
      t=0; %start time 
      %%%arrays to store results
j=1; %counter for arrays
t_array(1,t_end/t_sample)=zeros(1); t_array(1,j)=t; %time array and initial value
s_array(1,t_end/t_sample)=zeros(1); s_array(1,j)=s; %substrate array and initial value
p_array(1,t_end/t_sample)=zeros(1); p_array(1,j)=p; %product array and initial value

     while t < t_end
       %calculate rxn propensities
    a = [k1*s*e k2*c k3*c];
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
      s=s-1;
      e=e-1;
      c=c+1;
      elseif(mu == 2)
      c=c-1;
      s=s+1;
      e=e+1;
      elseif(mu == 3) 
      p=p+1;
      c=c-1;
      e=e+1;
      end 
      %store/output time and species
    if t >= j*t_sample
        j=j+1;
        t_array(1,j)=j;
        s_array(1,j)=s;
        p_array(1,j)=p;
      sub(:,NC)=s_array(1,j);
      pro(:,NC)=p_array(1,j);
      time(:,NC)=t_array(1,j);
    end 
    end 
 %update waitbar
    if t >= k*alpha*t_sample
        k=k+1;
        waitbar(t/t_end)
    end     
end
close(w)
toc
figure(1) 
plot(t_array,p_array)
hold on
plot(t_array,s_array,'r')
        title('protein-substrate number averaged over cells')
        legend('[P]','[S]','Location','best');
        xlabel('time in seconds')
        ylabel('concentration S y P')
        %axis([0 0.5 0 10]);
