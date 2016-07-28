clear all
clc

A = [2 -1 0; -1 2.5 -1.5; 0 -1.5 3];

%Guesses for phi1 and phi2
v1 = (1/sqrt(2))*[1;0;1];
v2 = [0; 1; 0];

V = [v1, v2];


Ap = V'*A*V;

[a, lam] = eig(Ap);

S = V*a

lam

% for i=1:10
%    
%    Xi = V'*A*V
%    
%    [a, lam] = eig(Xi)
%    
%    V = V*a
%    
% end