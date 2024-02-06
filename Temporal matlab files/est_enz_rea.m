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

clear all

%set parameter value
k1=30;
k2=1;
k3=10;
%set ensemble size
N=3;
%time units
T=10;
%set initial condition for molecule count and for simulation time, for each member of ensemble
sum(1)=0;
sum(2)=0;
sum(3)=0;
sum(4)=0;
su=5;en=1;co=0;pr=0;
S=zeros(T,N);
t=zeros(T,N);
S(1,:)=5;
t(1,:)=0;
E=zeros(T,N);
%Et=zeros(T,N);
E(1,:)=1;
%Et(1,:)=0;
C=zeros(T,N);
%Ct=zeros(T,N);
%C(1,:)=0;
%Ct(1,:)=0;
P=zeros(T,N);
%Pt=zeros(T,N);
%P(1,:)=0;
%Pt(1,:)=0;
%generate ensemble
for j=1:N
%run ten steps of the simulation algorithm
for i=1:T
    
%calculate propensity
%if C(i,j)>0
      a1=k1*su*en;
      a2=k2*co;
      a3=k3*co;
      a0=a1+a2+a3;
      %r1=rand(1,1);
      r2=rand(1);
      %tau=-(1./a0)*log(r1);
      yr2=r2*a0;    
%end
sum(2)=sum(1)+a1;
sum(3)=sum(2)+a2;
sum(4)=sum(3)+a3;
  for k=2:4
      if (sum(k) >= yr2) && (sum(k-1) < yr2) 
          mu=k-1;
      end 
  end
      if (mu == 1)
      su=su-1;
      en=en-1;
      co=co+1;
      S(i,j)=max(su,0);
      E(i,j)=max(en,0);
      C(i,j)=max(co,0);
      end
      if (mu == 2)
      co=co-1;
      su=su+1;
      en=en+1;
      S(i,j)=max(su,0);
      E(i,j)=max(en,0);
      C(i,j)=max(co,0);    
       end  
      if (mu == 3) 
      pr=pr+1;
      co=co-1;
      en=en+1;
      P(i,j)=max(pr,0);
      E(i,j)=max(en,0);
      C(i,j)=max(co,0);    
      end  
%update counter
%i=i+1;
%update time
t(i,j)=t(i,j)+(-(1/a0)*log(rand(1)));
%decrement state
%X1(i,j)=max(X1(i-1,j)-1,0);


end


end

%collect vector of reaction time events 
rxntimes=[];
for j=1:N
rxntimes=union(rxntimes, t(:,j));
end

%interpolate molecule counts for plotting (otherwise the ensemble of staircase graphs lie
%on top of one another)
for j=1:N
Sinterp(:,j)=interp1(t(:,j), S(:,j), rxntimes);
end
for j=1:N
Einterp(:,j)=interp1(t(:,j), E(:,j), rxntimes);
end
for j=1:N
Cinterp(:,j)=interp1(t(:,j), C(:,j), rxntimes);
end
for j=1:N
Pinterp(:,j)=interp1(t(:,j), P(:,j), rxntimes);
end

Sinterp(isnan(Sinterp)) = 0;
Smean=mean(Sinterp');
Einterp(isnan(Einterp)) = 0;
Emean=mean(Einterp');
Cinterp(isnan(Cinterp)) = 0;
Cmean=mean(Cinterp');
Pinterp(isnan(Pinterp)) = 0;
Pmean=mean(Pinterp');

%produce Figure 7.45A
figure(1)
%subplot(1,3,1)
colormap(gray)
c=[0 0 0;0.4 0.4 0.4;0.8 0.8 0.8];
set(0,'DefaultAxesColorOrder',c)
hold on
plot(rxntimes, Smean, 'k', 'Linewidth', 3)
hold on
plot(rxntimes, Emean, 'k', 'Linewidth', 3)
hold on
plot(rxntimes, Cmean, 'k', 'Linewidth', 3)
hold on
plot(rxntimes, Pmean, 'k', 'Linewidth', 3)
hold on
for j=1:N
plot(t(:,j), S(:,j), 'b', [0.5 0.5 0.5])
end
hold on
for j=1:N
plot(t(:,j), E(:,j), 'r', [0.5 0.5 0.5])
end
hold on
for j=1:N
plot(t(:,j), C(:,j), 'g', [0.5 0.5 0.5])
end
hold on
for j=1:N
plot(t(:,j), P(:,j), 'm', [0.5 0.5 0.5])
end
hold on
set(gca,'fontsize',12)

axis([0 5 0 10])
ylabel('Number of Molecules')
xlabel('Time')

%repeat for larger ensemble sizes

%set ensemble size
% N=10
% X2=zeros(10,N);
% t2=zeros(10,N);
% X2(1,:)=10;
% t2(1,:)=0;
% 
% for j=1:N
% 
% for i=1:10
%     
% %calculate propensity
% if X2(i,j)>0
% a=k*X2(i,j);
% end
% 
% %update counter
% i=i+1;
% 
% %update time
% t2(i,j)=t2(i-1,j)+log(1/rand(1))/a;
% 
% %decrement state
% X2(i,j)=max(X2(i-1,j)-1,0);
% 
% 
% end
% 
% 
% end
% 
% rxntimes=[];
% for j=1:N
% rxntimes=union(rxntimes, t2(:,j));
% end
% 
% 
% for j=1:N
% X2interp(:,j)=interp1(t2(:,j), X2(:,j), rxntimes);
% end
% 
% X2interp(isnan(X2interp)) = 0;
% X2mean=mean(X2interp');
% 
% %produce Figure 7.45B
% figure(1)
% hold on
% subplot(1,3,2)
% hold on
% for j=1:N
% plot(t2(:,j), X2(:,j), 'color', [0.5 0.5 0.5])
% end
% plot(rxntimes, X2mean, 'k','Linewidth', 3)
% axis([0 5 0 10])
% set(gca,'fontsize',12)
% ylabel('Number of Molecules')
% xlabel('Time')
% 
% 
% %set ensemble size
% N=30
% X3=zeros(10,N);
% t3=zeros(10,N);
% X3(1,:)=10;
% t3(1,:)=0;
% 
% for j=1:N
% 
% for i=1:10
%     
% %calculate propensity
% if X3(i,j)>0
% a=k*X3(i,j);
% end
% 
% %update counter
% i=i+1;
% 
% %update time
% t3(i,j)=t3(i-1,j)+log(1/rand(1))/a;
% 
% %decrement state
% X3(i,j)=max(X3(i-1,j)-1,0);
% 
% 
% end
% 
% 
% end
% 
% rxntimes=[];
% for j=1:N
% rxntimes=union(rxntimes, t3(:,j));
% end
% 
% 
% for j=1:N
% X3interp(:,j)=interp1(t3(:,j), X3(:,j), rxntimes);
% end
% 
% X3interp(isnan(X3interp)) = 0;
% X3mean=mean(X3interp');
% 
% %produce Figure 7.45C
% figure(1)
% hold on
% subplot(1,3,3)
% hold on
% for j=1:N
% plot(t3(:,j), X3(:,j), 'color', [0.5 0.5 0.5])
% end
% plot(rxntimes, X3mean, 'k','Linewidth', 3)
% axis([0 5 0 10])
% set(gca,'fontsize',12)
% ylabel('Number of Molecules')
% xlabel('Time')