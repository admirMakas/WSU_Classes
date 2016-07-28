clear all
clc

K=[1 -1 0;
    -1 2 -1;
    0 -1 1];
M=eye(3);

%K matrix is singular so needs to be shifted. In this case the diagonal
%will be shifted by 1.

Kbar = K+M;

%Need to guess X1bar
x1bar = [1 0; 1 1; 1 0];

for i=1:10
    
    %Get X2bar
    x2bar = (Kbar\M)*x1bar
    
    %get reduced K
    k1 = x2bar'*Kbar*x2bar
    
    %get reduced M
    m1 = x2bar'*M*x2bar
    
    [q, lam] = eig(m1\k1)
    
    %q = q(:, [2,1])
    
    lam=diag(sort(diag(lam)))
    
    x1bar = x2bar*q
    %x1bar = x1bar/[norm(x1bar(:,1)) 0; 0 norm(x1bar(:,2))]

end