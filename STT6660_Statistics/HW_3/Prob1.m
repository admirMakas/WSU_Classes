clear all
clc

% x = [27 45 41 19 35 39];
% y = [57 64 80 46 62 72];

x = [1 2 3 4 5];
y = [1 1 2 2 4];

%%
% part a - make scatter plot of data

%scatter(x,y)

%%
% part b - find b1 and b0, write down equation
avg_x = mean(x)
avg_y = mean(y);

n = size(x,2);

%scatter(x,y)

Sxy = sum(x.*y) - (sum(x)*sum(y))/n
Sxx = sum(x.*x) - (sum(x)^2)/n

b1 = Sxy/Sxx

b0 = avg_y - b1*avg_x

%%
% part c - make scatter plot of residuals vs predictor

y_hat = b1.*x + b0

e = y-y_hat

% scatter(x, e)

%%
% part d - make scatter plot of residuals vs fitted values. Compute SSE

%scatter(y_hat, e)

SSE = sum(e.^2)

%%
% part h - find 95% CI for B1

prob = .975;

Zt = tinv(prob, 3)

Probt = tcdf(1.237, 7)