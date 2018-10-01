clear;
%initial condition
p=1;
x0=3.1;
m=1;
L=25;

N=30;
%initial condition of gradient method with alpha=1/L
xG1=ones(1,2)*x0;
%initial condition of gradient method with alpha=2/(m+L)
xG2=ones(1,2)*x0;
%initial condition of Nesterov's accelerated method
xN=ones(1,2)*x0;
%initial condition of Heavy-ball method
xH=ones(1,2)*x0;
for k=2:N;
    %Gradient method with alpha=1/L
     if xG1(k)<1
         g = 25*xG1(k);
      elseif xG1(k)<2 
         g = xG1(k)+24;
      else
         g = 25*xG1(k)-24;
     end
    xG1(k+1)=xG1(k)-(1/L)*g;
    %Gradient method with alpha=2/(m+L)
      if xG2(k)<1
         g = 25*xG2(k);
      elseif xG2(k)<2 
         g = xG2(k)+24;
      else
         g = 25*xG2(k)-24;
     end
    xG2(k+1)=xG2(k)-(2/(m+L))*g;
    %Nesterov's method
    bN=(sqrt(L)-sqrt(m))/(sqrt(L)+sqrt(m));
     if (1+bN)*xN(k)-bN*xN(k-1)<1
         g = 25*((1+bN)*xN(k)-bN*xN(k-1));
      elseif  (1+bN)*xN(k)-bN*xN(k-1)<2
         g = (1+bN)*xN(k)-bN*xN(k-1)+24;
      else
         g = 25*((1+bN)*xN(k)-bN*xN(k-1))-24;
     end
    xN(k+1)=(1+bN)*xN(k)-bN*xN(k-1)-(1/L)*g;
    %Heavy-ball method
    aH=4/(sqrt(L)+sqrt(m))^2;
    bH=bN^2;
     if xH(k)<1
         g = 25*xH(k);
      elseif xH(k)<2 
         g = xH(k)+24;
      else
         g = 25*xH(k)-24;
     end
    xH(k+1)=(1+bH)*xH(k)-bH*xH(k-1)-aH*g;
end

%plot results
plot(0:N-1, abs(xG1(2:N+1)),'b*', 0:N-1, abs(xG2(2:N+1)),'g--', 0:N-1, abs(xN(2:N+1)),'r-.', 0:N-1, abs(xH(2:N+1)),'k')
legend('Gradient Method with \alpha=1/L', 'Gradient Method with \alpha=2/(m+L)', 'Nesterov Method', 'Heavy-ball Method')
title(['L=',num2str(L),', m=',num2str(m)])
xlabel('steps(k)');
ylabel('|x(k)|');







