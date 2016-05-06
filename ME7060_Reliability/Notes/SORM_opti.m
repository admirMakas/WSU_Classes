function HL_opti


% run optimization within the U space



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FORM to find MPP
u0=[0 0 ];
lb=[-5 -5 ]
ub=[5 5]

options = optimset('display','iter');
[u,fval]=fmincon(@objf,u0,[],[],[],[],lb,ub,@constf,options)
beta=norm(u);
u;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SORM 

% B matrix

[guk,dguk,d2guk]=limit_state_gu(u);
B=diag(d2guk)./norm(dguk)

% H matrix

H0=[-dguk'./norm(dguk);0 1]'
[m,n]=size(H0)
H=zeros(m,n)
R=zeros(n,n)
for k=1:n
    ri=H0(:,k)
    for i=1:k-1
        R(i,k)=H(:,i)'*H0(:,k)
        ri=ri-R(i,k)*H(:,i)
    end
    R(k,k)=norm(ri)
    H(:,k)=ri/R(k,k)
end
H=[H(2:n,:);H(1,:)]


% Curvature k
K=H*B*H'

k1=K(1,1)

% Compute pf  (Breitung's formula)

pf=normcdf(-beta)*(1+k1*beta)^(-1/2)



% Compute pf  (Tvedt's formula)

A1=normcdf(-beta)*(1+k1*beta)^(-1/2)
A2=(beta*normcdf(-beta)-normpdf(beta))*((1+beta*k1)^(-1/2)-(1+(beta+1)*k1)^(-1/2))
A3=(beta+1)*(beta*normcdf(-beta)-normpdf(beta))*((1+beta*k1)^(-1/2)-real((1+(beta+1)*k1)^(-1/2)))
pf=A1+A2+A3



function y=objf(u)

y=norm(u);


function [c, ceq]=constf(u)

c=[];
ceq=limit_state_gu(u);
