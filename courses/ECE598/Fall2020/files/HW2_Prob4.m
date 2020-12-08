%%
clear;
% Problem Setup

A = [0.99 0.01 0;0.01 0.95 0.01;0.5 0.14 0.97];
B = [1 0.1;0 0.2;1.1 0.7];
Q = blkdiag(1, 1, 0.5);
R = blkdiag(1.3, 0.7);
W = blkdiag(0.3, 0.15, 0.1);
gamma = 0.98;


%%
%Generate data for LSPI
%Notice the data generated for LSPI does not require knowing controller K

%Specify the sample size
N=200000;

%Generate x and u from uniform distribution [-C, C]
%It is possible to use other behavior policy to sample u_data
%the uniform distribution does provide enough exploration
C=500;
x_data=rand(3,N)*2*C-C;
u_data=rand(2,N)*2*C-C;
% generate x1=Ax+Bu+w
for i=1:N
    x1_data(:,i)=A*x_data(:,i)+B*u_data(:,i)+mvnrnd([0 0 0],W)';
end







%%
%Problem 4

%solve optimal K from optimal Bellman equation (DARE) 
%K_opt will be used for comparison
P_opt=dare(sqrt(gamma)*A,sqrt(gamma)*B,Q,R);
K_opt=gamma*inv(R+gamma*B'*P_opt*B)*B'*P_opt*A;

%solve optimal K using approximate PI with LSTDQ for Q estimation

Num_iter=10;

%Initialize K
%One needs a stabilizing controller here
K=[0 0.1 0; 0.1 0 0.1];

for j=1:Num_iter

    A1=zeros(16);
b1=zeros(16,1);
for i=1:N
    % read the i-th data point
    x=x_data(:,i);
    u=u_data(:,i);
    %calculate the cost
    c=x'*Q*x+u'*R*u;
    % read x1 from x1_data
    xN=x1_data(:,i);
% generate uN=-K*xN
    uN=-K*xN;
    
XN=[xN;uN]*([xN;uN]');
X=[x;u]*([x;u]');
%calculate the feature phi(x1,-K*x1)
xxN=XN(triu(ones(size(XN)))~=0);
%calculate the feature phi(x,u)
xx=X(triu(ones(size(X)))~=0);
%augment the feature with the residue term
xx=[xx;1];
xxN=[xxN;1];

A1=A1+xx*((xx-gamma*xxN)');
b1=b1+xx*c;
end
%use LSTDQ to estimate the weight vector
O_LSTDQ=pinv(A1)*(b1);
    
%rewrite the weight into the Q matrix
%read Q22 from the weight
Q22=[O_LSTDQ(10) O_LSTDQ(14)/2;O_LSTDQ(14)/2 O_LSTDQ(15)];
%read Q12 from the weight
Q12=[O_LSTDQ(7:9) O_LSTDQ(11:13)]/2;
%update K matrix

K=inv(Q22)*Q12';

end

%compare K with optimal K
(K-K_opt)./K_opt
K
K_opt

