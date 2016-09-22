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

figure()
scatter(x, y, 90, 'fillcolor', 'k')
axis([-1 6 0 5])
title('Scatter plot of X vs Y', 'fontsize', 16)
xlabel('Amount of Drug %', 'fontsize', 16)
ylabel('Time (sec)', 'fontsize', 16)

Sxy = sum(x.*y) - (sum(x)*sum(y))/n
Sxx = sum(x.*x) - (sum(x)^2)/n

b1 = Sxy/Sxx

b0 = avg_y - b1*avg_x

%%
% part c - make scatter plot of residuals vs predictor

y_hat = b1.*x + b0

e = y-y_hat

figure()
scatter(x, e, 90, 'fillcolor', 'k')
axis([-1 6 -2 2])
title('Error vs Predictor', 'fontsize', 16)
xlabel('Perdictor, Amount of Drug %', 'fontsize', 16)
ylabel('Error', 'fontsize', 16)

%%
% part d - make scatter plot of residuals vs fitted values. Compute SSE

figure()
scatter(y, e, 90, 'fillcolor', 'k')
axis([0 5 -1.5 1.5])
title('Error vs Fitted Values', 'fontsize', 16)
xlabel('Fitted Values, Time (sec)', 'fontsize', 16)
ylabel('Error', 'fontsize', 16)

SSE = sum(e.^2)

%%
% part h - find 95% CI for B1

prob = .975;

Zt = tinv(prob, 3)

Probt = tcdf(1.237, 7)