clear;
%specify the problem dimension
p=10;
%you should specifiy  the value of L by yourself
L=1000;
%you should specify the value of m by yourself
m=0.139;


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