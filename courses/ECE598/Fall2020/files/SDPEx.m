clear;
% set A matrix, change this based on your need
A=eye(3)*0.9;
%implement the LMI
cvx_begin sdp quiet
    variable P(3,3) semidefinite;     
     minimize(0)
    subject to
    A'*P*A-P<=-0.001*eye(3); 
 P>=0.001*eye(3)
cvx_end

P

%another way to break homogeneity
cvx_begin sdp quiet
    variable P(3,3) semidefinite;     
     minimize(0)
    subject to
    A'*P*A-P<=0; 
 P>=0;
 trace(P)==1;
cvx_end

P