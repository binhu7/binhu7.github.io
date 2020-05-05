%%
clear;
% Problem Setup
A=[0.99 0.01 0; 0.01 0.98 0.01; 0.5 0.12 0.97];
B=[1 0.1; 0 0.1; 0 0.1];
Q=eye(3);
R=eye(2);
W=blkdiag(0.1,0.05,0.1);
gamma=0.98;

%%
%Problem 5(a)
%Specify the controller K
K=[0 0.1 0; 0.1 0 0.1];
%Calculate the state value function 
P=dlyap(sqrt(gamma)*(A-B*K)',Q+K'*R*K)
%Calculate Q from P
Op=gamma*[A B]'*P*[A B]+blkdiag(Q,R);
%Directly calculate Q from the Bellman equation for Q
Op1=dlyap(sqrt(gamma)*[A B;-K*A -K*B]',blkdiag(Q,R))
%Compare the above two approaches
Op-Op1

%%
%Generate data for LSTDQ
%Notice the data generated for LSTDQ does not require knowing controller K

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
%Problem 5(b)
%LSTD for Q
%LSTD for state-action value function
%Generate the weight of the true Q for comparison
Opt=Op*2;
for i=1:5
Opt(i,i)=Opt(i,i)/2;
end
%Copy the quadratic part into the weight
Ostar=Opt(triu(ones(size(Opt)))~=0);




%Specify the controller K
K=[0 0.1 0; 0.1 0 0.1];

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
O_LSTDQ=inv(A1)*(b1);
%compare with the true weight
(O_LSTDQ(1:15)-Ostar)./Ostar



%%
%Problem 5(c)

%solve optimal K from optimal Bellman equation (DARE) 
%K_opt will be used for comparison
P_opt=dare(sqrt(gamma)*A,sqrt(gamma)*B,Q,R);
K_opt=gamma*inv(R+gamma*B'*P_opt*B)*B'*P_opt*A;

%solve optimal K using approximate PI with LSTDQ for Q estimation

Num_iter=10;

%Initialize K
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

