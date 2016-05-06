function HL_opti


% run optimization within the U space




u0=[0 0 0];
lb=[-5 -5 -5]
ub=[5 5 5]

options = optimset('display','iter');
[u,fval]=fmincon(@objf,u0,[],[],[],[],lb,ub,@constf,options)

function y=objf(u)

y=norm(u);


function [c, ceq]=constf(u)

c=[];
ceq=limit_state_gu(u);
