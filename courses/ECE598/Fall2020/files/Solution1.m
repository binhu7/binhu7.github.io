%%
%HW1
%%%%%%%%%%%%%%%%%
%Problem 1(a)
%Set up the A matrix
A=[-0.5 -0.85 0.33;
    -0.23 0.71 0.79;
    0.41 -0.1 -0.21];
%implement the LMI
cvx_begin sdp quiet
    variable P(3,3) semidefinite;     
     minimize(0)
    subject to
    % break homegeneity by choosing epsilon to be 0.001
    A'*P*A-P<=-0.001*eye(3); 
 P>=0.001*eye(3)
cvx_end

P

clear P;
cvx_begin sdp quiet
    variable P(3,3) semidefinite;     
     minimize(0)
    subject to
    % break homegeneity by choosing epsilon to be 0.001
    A'*P*A-P<=0; 
 P>=0;
 trace(P)==1;
cvx_end

P


clear P;

%%%%%%%%%%%%%%%%%%%%%%%
%%
% Problem 1(b)
eig(A);
%calculate the spectral radius of A
rho1=max(abs(ans))

%A=eye(3)*0.9;
%implement the LMI
cvx_begin sdp quiet
    variable P(3,3) semidefinite;     
     minimize(0)
    subject to
    A'*P*A-0.97293^2*P<=-0.001*eye(3); 
 P>=0.001*eye(3)
cvx_end

P
%%
% Problem 1(c)
clear;
% specify the value of (m,L)here
m=1;
L=100;
%specify the value of alpha here
% alpha has to be smaller than 2/L
alpha=1/L;
rho=max(abs(1-m*alpha),abs(1-L*alpha));
cvx_begin sdp quiet
    variables lambda r2; 
    %r2 is the value of rho^2 in the LMI and can be directly minimized
     minimize(r2)
    subject to
   [1-r2 -alpha;-alpha alpha^2]<=lambda*[2*m*L -(m+L);-(m+L) 2]; 
 lambda>=0;
 %r2 has to be greater than 0 and smaller than 1
 r2>=0;
 r2<=1;
cvx_end
%the value of lambda
lambda
%compare r2 with the theoretical value rho^2
% you can find these two values are quite close
r2-rho^2

%%
%Problem 2(c)
% we can get rid of the Kronecker product in the implementation
% specify the value of (m,L)here
clear;
m=1;
L=100;
%specify the value of n
n=20;
%specify the value of alpha here
% specify the value of alpha
alpha=1/(3*L);
% specify the value of rho
r2=1-min(1/(3*n),m/(3*L));
%specify e_i as I(:,i)
I=eye(n);
%specify e
e=ones(n,1);


%specify A_i and B_i
for i=1:n
ei=I(:,i);
    A(:,:,i)=[eye(n)-ei*ei' zeros(n,1); -(alpha/n)*(e-n*ei)' 1];
B(:,:,i)=[ei*ei';-alpha*ei'];
end
%specify the value of C
C=[zeros(1,n) 1];
%specify supply rates
%use Assumption 2 on Page 14 of Lecture Note 4
M=[1 zeros(1,n); 0 e'/n]'*[2*m*L -(m+L);-(m+L) 2]*[1 zeros(1,n); 0 e'/n];
X0=[C zeros(1,20);zeros(20,21) eye(20)]'*M*[C zeros(1,20);zeros(20,21) eye(20)];
%use Assumption 1 on page 14 of Lecture Note 4
for i=1:n
    ei=I(:,i);
M=[1 zeros(1,n); 0 ei']'*[2*m*L -(m+L);-(m+L) 2]*[1 zeros(1,n); 0 ei'];
X(:,:,i)=[C zeros(1,20);zeros(20,21) eye(20)]'*M*[C zeros(1,20);zeros(20,21) eye(20)];
end
%choice 1
%formulate the LMI
cvx_begin sdp quiet
   variable P(n+1,n+1) semidefinite;
    variable lambda(n) nonnegative;
    variable lambda0 nonnegative;
    X_sum=lambda0*X0;
    L_sum=zeros(41,41);
    for i=1:n
    X_sum=X_sum+lambda(i)*X(:,:,i);
    L_sum= L_sum+[A(:,:,i)'*P*A(:,:,i)-r2*P A(:,:,i)'*P*B(:,:,i);B(:,:,i)'*P*A(:,:,i) B(:,:,i)'*P*B(:,:,i)];
    end
    L_sum=L_sum/n;
     minimize(0)
    subject to
    L_sum<=X_sum;
%break homogeneity
    trace(P)==1;
    %p2>=0;
cvx_end
P

%choice 2
%Enforce P to be diagonal in  the LMI and reduce the number of lambda
cvx_begin sdp quiet
   % variable P(n+1,n+1) semidefinite;
   variables p1 p2;
   P=blkdiag(p1*eye(n),p2);
    %variable lambda(n) nonnegative;
   % variable lambda0;
   variables lambda0 lambda1;
    %X_sum=lambda0*X0;
    L_sum=zeros(41,41);
    X_sum=zeros(41,41);
    for i=1:n
    X_sum=X_sum+X(:,:,i);
    L_sum= L_sum+[A(:,:,i)'*P*A(:,:,i)-r2*P A(:,:,i)'*P*B(:,:,i);B(:,:,i)'*P*A(:,:,i) B(:,:,i)'*P*B(:,:,i)];
    end
    L_sum=L_sum/n;
     minimize(0)
    subject to
    L_sum<=lambda1*X_sum+lambda0*X0;
    lambda0==0;
    lambda1>=0;
%break homogeneity
    %p1>=0.001;
    p1==1/L;
    p2>=0;
cvx_end
P



%choice 3
%Enforce p1=2/(3*L), p2=1/alpha, lambda1=1/(L*n)
cvx_begin sdp quiet
   % variable P(n+1,n+1) semidefinite;
   variables p1 p2;
   P=blkdiag(p1*eye(n),p2);
    %variable lambda(n) nonnegative;
   % variable lambda0;
   variables lambda0 lambda1;
    %X_sum=lambda0*X0;
    L_sum=zeros(41,41);
    X_sum=zeros(41,41);
    for i=1:n
    X_sum=X_sum+X(:,:,i);
    L_sum= L_sum+[A(:,:,i)'*P*A(:,:,i)-r2*P A(:,:,i)'*P*B(:,:,i);B(:,:,i)'*P*A(:,:,i) B(:,:,i)'*P*B(:,:,i)];
    end
    L_sum=L_sum/n;
     minimize(0)
    subject to
    L_sum<=lambda1*X_sum+lambda0*X0;
    lambda0==0;
    lambda1==1/L/n;
%break homogeneity
    %p1>=0.001;
    p1==2/(3*L);
    p2==1/alpha;
    %p2>=0;
cvx_end
P






%%
% Problem 3(c)
% we can get rid of the Kronecker product in the implementation
% specify the value of (m,L)here
clear;
m=1;
L=100;
%specify the value of alpha here
% specify the value of alpha
alpha=1/L;
%specify the value of beta
beta=(sqrt(L)-sqrt(m))/(sqrt(L)+sqrt(m));
% specify the value of rho
r2=1-sqrt(m/L);

%specify A and B
A=[1+beta -beta;1 0];
B=[-alpha;0];
%specify X1, X2, and X
X1=0.5*[beta^2*m -beta^2*m -beta;-beta^2*m beta^2*m beta;-beta beta alpha*(2-L*alpha)];
X2=0.5*[(1+beta)^2*m -beta*(1+beta)*m -(1+beta);-beta*(1+beta)*m beta^2*m beta;-(1+beta) beta alpha*(2-L*alpha)];
X=r2*X1+(1-r2)*X2;
%formulate the LMI
cvx_begin sdp quiet
    variable P(2,2) semidefinite;
     minimize(0)
    subject to
    [A'*P*A-r2*P A'*P*B;B'*P*A B'*P*B]<=X;
cvx_end
P
[A'*P*A-r2*P A'*P*B;B'*P*A B'*P*B]-X




