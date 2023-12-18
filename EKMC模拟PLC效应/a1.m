 clc
clear
fnameetotal=sprintf('etotal_%d .vtk',1);
    outetotal=fopen(fnameetotal,'w');
    fnamebtotal=sprintf('btotal_%d .vtk',1);
    outbtotal=fopen(fnamebtotal,'w');
%³õÊ¼Ìõ¼þ
ep0(32)=0;
ep0(:)=0.002;
ep0v(32)=2e-2;
ta0(32)=0;

for i=1:32
    
    b0(i)=real(5.2*log(ep0v(i)/0.94e-4+sqrt((ep0v(i)/0.94e-4).^2+1))+62.5+...
        180*(1-exp(-15.8*ep0(i)))*(1+0.42*log(ep0(i)).^0.57)+...
        6.44*1.38*(1-exp(-0.4*(ep0(i)).^0.67*ta0(i).^0.33)));
   
end

x=-15;
dt=1e-4;
for i=1:32 
   b1(i)=b0(i)+x;
   x=x+1;
end

%%%%%%%%%%%
step=1;
t0=0;
dt=1e-4;
btotaL1=0;
btotaL0=170;
while (step<=3000)
b1
ep0
    for i=1:32
    syms k
   ep0v(i)=real(vpa(solve(b1(i)==5.2*log(k/(0.94e-4)+sqrt((k/0.94e-4).^2+1))+62.5+...
        180*(1-exp(-15.8*ep0(i)))*(1+0.42*log((ep0(i)).^0.57))+...
        6.44*1.38*(1-exp(-0.4*((ep0(i)).^0.67)*(ta0(i).^0.33))),k)));  
   ta0v(i)=1-ta0(i)/5e-4*ep0v(i);
   end


t1=t0+dt;

syms xx
for i=1:32
   ep1(i)=ep0(i)+ep0v(i)*dt;%2

   ta1(i)=ta0(i)+ta0v(i)*dt;%3

   btotalv(i)=double(((2e-2)-1/32*int(ep0v(i),xx,0,32))*7e4);%4

   m(i)=double(btotaL0+btotalv(i)*dt);%5
   btotaL1=btotaL1+m(i);
end
btotaL1=1/32*btotaL1;

%6
F=double(12*btotaL1/(1+2e-2*t1));
%7
F
i=1;
    for x=-15:16
   A(i)=12e-4/(1+2e-2*t1)*x+12/(1+2e-2*t1);
   i=i+1;
    end
etotaL1=2e-2*t1;
for i=1:32
     b2(i)=real(F*(1+ep1(i)+b1(i)/7e4)/A(i));
end



        fprintf(outetotal,'%14.6e\n',double(etotaL1));
        fprintf(outbtotal,'%14.6e\n',double(btotaL1));
 for i=1:32
     ep0(i)=ep1(i);
     ta0(i)=ta1(i);
     b1(i)=b2(i);
 end

 t0=t1;
 btotaL0=btotaL1;
% 
%     for i=1:32
%     syms kk
%    kk=(solve(b1(i)==5.2*log(kk/0.94e-4+sqrt((kk/0.94e-4).^2+1))+62.5+...
%         180*(1-exp(-15.8*ep0(i)))*(1+0.42*log(ep0(i)).^0.57)+...
%         6.44*1.38*(1-exp(-0.4*(ep0(i)).^0.67*ta0(i).^0.33)),kk));
%     kk
%     vpa(kk)
%     double(kk)
%     
%      
%     ta0v(i)=1-ta0(i)/5e-4*ep0v(i);
%    end
fprintf('done step: %5d\n',step);
 step=step+1;
end
        
 fclose(outetotal);
 fclose(outbtotal);