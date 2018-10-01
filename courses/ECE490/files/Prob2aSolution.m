clear;
%specify the problem dimension
p=500;
%you should specifiy  the value of L by yourself
L=10;
%you should specify the value of m by yourself
m=1;


%compute the condition number
kappa=L/m;

%generate Q
Q=randn(p,p); [A, R]=qr(Q);
D=rand(p,1);D=10.^D;Dmin=min(D);Dmax=max(D);
D=(D-Dmin)/(Dmax-Dmin);
D=m+D*(L-m);
Q=A'*diag(D)*A;

%generate vector q
q=randn(p,1);

%generate scalar r
r=randn(1,1);


N=350;
%optimal solution x*
xstar=-inv(Q)*q;
fstar=(0.5)*xstar'*Q*xstar+q'*xstar+r;
%initial condition of gradient method with alpha=1/L
xG1=ones(p,2);
fG1=(0.5)*xG1(:,2)'*Q*xG1(:,2)+q'*xG1(:,2)+r;
%initial condition of gradient method with alpha=2/(m+L)
xG2=ones(p,2);
fG2=(0.5)*xG2(:,2)'*Q*xG2(:,2)+q'*xG2(:,2)+r;
%initial condition of Nesterov's accelerated method
xN=ones(p,2);
fN=(0.5)*xN(:,2)'*Q*xN(:,2)+q'*xN(:,2)+r;
%initial condition of Heavy-ball method
xH=ones(p,2);
fH=(0.5)*xH(:,2)'*Q*xH(:,2)+q'*xH(:,2)+r;
for k=2:N;
    %Gradient method with alpha=1/L
    xG1(:,k+1)=xG1(:,k)-(1/L)*(Q*xG1(:,k)+q);
    fG1(k)=(0.5)*xG1(:,k+1)'*Q*xG1(:,k+1)+q'*xG1(:,k+1)+r;
    %Gradient method with alpha=2/(m+L)
    xG2(:,k+1)=xG2(:,k)-(2/(m+L))*(Q*xG2(:,k)+q);
    fG2(k)=(0.5)*xG2(:,k+1)'*Q*xG2(:,k+1)+q'*xG2(:,k+1)+r;
    %Nesterov's method
    bN=(sqrt(L)-sqrt(m))/(sqrt(L)+sqrt(m));
    xN(:,k+1)=(1+bN)*xN(:,k)-bN*xN(:,k-1)-(1/L)*(Q*((1+bN)*xN(:,k)-bN*xN(:,k-1))+q);
    fN(k)=(0.5)*xN(:,k+1)'*Q*xN(:,k+1)+q'*xN(:,k+1)+r;
    %Heavy-ball method
    aH=4/(sqrt(L)+sqrt(m))^2;
    bH=bN^2;
    xH(:,k+1)=(1+bH)*xH(:,k)-bH*xH(:,k-1)-aH*(Q*xH(:,k)+q);
    fH(k)=(0.5)*xH(:,k+1)'*Q*xH(:,k+1)+q'*xH(:,k+1)+r;
end
%plot results
semilogy(1:N, fG1-fstar,'b:', 1:N, fG2-fstar,'g--', 1:N, fN-fstar,'r-.', 1:N, fH-fstar,'k')
legend('Gradient Method with \alpha=1/L', 'Gradient Method with \alpha=2/(m+L)', 'Nesterov Method', 'Heavy-ball Method')
title(['p=',int2str(p),', L=',num2str(L),', m=',num2str(m)])
xlabel('steps(k)');
ylabel('log(f(x(k))-f(x*))');

