%function [muy, sigy, pf, X1] = MCS_Solver()
clear all
N=10000;

U1=rand(N,1);
U2=rand(N,1);
y=[];
for i=1:N
    u1=U1(i);
    u2=U2(i);
    %ui==Fu==Fx
    %Fx1=u1 and Fx2=u2
    
    %Uniform Distribution
    x1 = 10 + u1*(20-10); %variable P
    %Normal distribution
    x2 = 2 + 0.2*icdf('normal', u2, 0, 1); %variable w
    %Calculate Moment
    yi = 112.5*x2 + 7.5*x1;
    
    y=[y;yi];
end

muy = mean(y)
sigy = std(y)

nf=sum(y>350);
%Estimate failure probability
pf = nf/N;

%List memory addresses for the vector y components that are greater than
%350
ids = y>350;

X1 = 10 + U1*(20-10);
X2 = 2 + 0.2*icdf('normal', U2, 0, 1);
%figure
scatter(X1(ids), X2(ids), '.')
%scatter(X1(ids), x2(ids))
hold on
scatter(X1(~ids), X2(~ids), '*')
%end