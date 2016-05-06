function [muy,sigy,Pf]=MCS_solver(N)

N=30000

U1=rand(N,1);
U2=rand(N,1);
y=[];
for i=1:N
    u1=U1(i);
    u2=U2(i);
    %inverse transformation
    %ui==Fu==Fx
    %Fx1=u1  &  Fx2=u2
    
    
    % Uniform Distribution U[10, 20]
    x1=10+u1*(20-10); % variable P
    x2=2+0.2*icdf('normal',u2,0,1); %variable W~N(2.0,0.2)
    %x2=norminv(u2,2,0.2)
    yi=112.5*x2+7.5*x1; % Moment
    
        
    
    y=[y;yi];
end

muy=mean(y);
sigy=std(y);

nf=sum(y>350);
Pf=nf/N;


ids=y>350


X1=10+U1*(20-10);
x2=2+0.2*icdf('normal',U2,0,1);
close all
figure
scatter(X1(ids),x2(ids),'.')    
hold on
scatter(X1(~ids),x2(~ids),'.') 