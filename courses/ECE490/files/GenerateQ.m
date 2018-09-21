clear;
%specify the problem dimension
p=100;
%you should specifiy  the value of L by yourself
L=1000;
%you should specify the value of m by yourself
m=0.1;


%compute the condition number
kappa=L/m;

%generate Q
A=randn(p);
[V, D]=svd(A);
alpha=(kappa*D(p,p)/D(1,1))^(1/(p-1));
for i=1:p; a(i)=alpha^(p-i);end
D=D*diag(a);
D=D*(m/D(p,p));
Q=V'*D*V;

%generate vector q
q=randn(p,1);

%generate scalar r
r=randn(1,1);