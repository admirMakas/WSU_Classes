%%
% This code will supply the Z value given some probability

clear all
close all
clc

mu = 0;
sig = 1;

% mu = 170;
% sig = 30;
pd = makedist('Normal', mu, sig);

range = linspace(0,350,200);

cdf_val = cdf(pd, range);
pdf_val = pdf(pd, range);

% figure()
% plot(range, cdf_val)
% figure()
% plot(range, pdf_val)

r = iqr(pd);

prob = .9000;

Z = icdf(pd, prob)

%%
% This code will find the probability given some Z

Zp = [180 230];

Probi = cdf(pd, Zp)