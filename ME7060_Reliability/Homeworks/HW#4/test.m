clc
mux = 120;
sigx = 12;

xprime = 80.0402;

sigy = sqrt(log((sigx/mux)^2+1))
muy = log(mux)-0.5*sigy^2

% get pdf value for lognormal dist
pdf1 = pdf('lognorm', xprime, muy, sigy)

% get inverse CDF with normal distribution of the lognormal
% CDF value
cdf1 = icdf('norm', cdf('lognormal', xprime, muy, sigy), 0, 1)

% now can calculate equivalent normal mu and sig to be used.
sigxP = pdf('norm', cdf1, 0, 1)/pdf1

muxP = xprime - cdf1*sigxP