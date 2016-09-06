%%
%This code will supply the Z value given some probability
mu = 0;
sig = 1;
pd = makedist('Normal', mu, sig);

prob=.21;

Z=icdf(pd,prob)

%%
%This code will find the probability given some Z

Zi = [-1.647];

probi = cdf(pd,Zi);