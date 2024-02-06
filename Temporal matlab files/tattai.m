 %  MATLAB tattai.m  model from Thattai- van Oudenaarden 2001 paper.
%  SINGLE GENE MODEL
%       There are 4 reactions: r=mRNA,d=DNA,p=protein, in reactions =
%      stands for right arrow
% reaction d=r+d (c1)
% reaction r=0 (c2)
% reaction r=r+p (c3)
% reaction p=0 (c4)
% parameters of the reaction system, nens=nb of cells, timelimit= length
%    of running time, itime is length of file which gives time 
% dependence of proteins or averages
      itime=1000;
      nens=200;
      timelimit=30000;
      b=2;
%  parameters:reaction rates
      c1=0.01;
      c2=log(2)/120;
      c3=b*c2;
      c4=log(2)/3600;
      avp=0;
      avp2=0;
      tim=zeros(itime-1,1);
      npt=zeros(itime-1,1);
      si1=zeros(itime-1,1);
      si2=zeros(itime-1,1);
      si3=zeros(itime-1,1);
      npt2=zeros(itime-1,1);
      np1=zeros(itime-1,1);
      np12=zeros(itime-1,1);
        for  is=1:nens
%  intialization
      d=1;
      r=0;
      p=0;
      time=0;
      izt=1;
     while time < timelimit
      sum(1)=0;
      sum(2)=0;
      sum(3)=0;
      sum(4)=0;
      sum(5)=0;
      a1=c1*d;
      a2=c2*r;
      a3=c3*r;
      a4=c4*p;
      a0=a1+a2+a3+a4;
      r1=rand(1,1);
      r2=rand(1,1);
      tau=-(1./a0)*log(r1);
      yr2=r2*a0;
      sum(2)=sum(1)+a1;
      sum(3)=sum(2)+a2;
      sum(4)=sum(3)+a3;
      sum(5)=sum(4)+a4;
      for k=2:5
      if (sum(k) >= yr2) && (sum(k-1) < yr2) 
          mu=k-1;
      end 
        end
      if (mu == 1) 
      r=r+1;
      end
      if (mu == 2)
      r=r-1;
       end  
      if (mu == 3) 
      p=p+1;
       end  
      if (mu == 4) 
          p=p-1;
      end
       ztime=timelimit./itime;
       if time > izt*ztime
           time;
           izt;
           np1(izt)=np1(izt)+p;
           np12(izt)=np12(izt)+p^2;
             if (is ==23)
              si1(izt)=p;
                end
             if (is ==43)
              si2(izt)=p;
                end
             if (is ==63)
              si3(izt)=p;
                end
              
 %          np2(izt)=np2(izt)+r;
 %          np22(izt)=np22(izt)+r^2;
           izt=izt+1;
               end 
       time=time + tau;
          end
%   x is the vector containing steady state values for all cells
         np1;
        x(is)=p;
        xr(is)=r;
       avp=avp+p;
       avp2=avp2+p.^2;
         end
%   steady state measures for proteins averaged over all cells 
       avp=avp./nens
       avp2=avp2./nens;
       var=avp2-avp.^2;
       standarddev=sqrt(var)
       fano=var./avp
       cvar=sqrt(var./avp^2)
%  time behavior of proteins averaged over cells
       for jt=1:itime-1
       tim(jt)=ztime*jt;
        end
       np=np1./nens;
       np12=np12./nens;
       coeff=(np12-np1.^2);
       fao=coeff./np1;
       cvars=fao./np1;
%       coeff=sqrt(coeff);
       tim1=tim([1:itime-10],[1]);
       fao1=fao([1:itime-10],[1]);
       cvars1=cvars([1:itime-10],[1]);
       npt1=np([1:itime-10],[1]);
       si11=si1([1:itime-10],[1]);
       si21=si2([1:itime-10],[1]);
       si31=si3([1:itime-10],[1]);
       %subplot(3,2,1), 
       figure(1)
       plot(tim,np)
       title('protein number averaged over cells, b=2')
       xlabel('time in seconds')
       legend('200 cells',4)
        %subplot(3,2,2),
        figure(2)
        plot(tim1,npt1)
        title('protein number: average and 3 single cells')
        xlabel('time in seconds')
        hold on
        %subplot(3,2,2), 
        plot(tim1,si11,'r')
        hold on
        %subplot(3,2,2), 
        plot(tim1,si21,'g')
        hold on
        %subplot(3,2,2), 
        plot(tim1,si31,'m')
%        subplot(3,2,3),  plot(tim1,cvars1)
%        title('protein coeff. var. squared, averaged over cells, b=2')
%        xlabel('time in seconds')
%        legend('200 cells',4)
%        subplot(3,2,4),  plot(tim1,fao1)
%        title('protein Fano factor, averaged over cells, b=2')
%        xlabel('time in seconds')
%        legend('200 cells',4)
% %       hold off
% %  histogram of steady state protein numbers
%        subplot(3,2,5)
figure(3)
         hist(x,10)
        title('histogram of steady state protein number')
%         hold off
%        subplot(3,2,6)
figure(4)
        hist(xr,50)
        title('histogram of steady state mRNA number')  
%         hold off

%
