clear all
clc

x = [1 0 2 0 3 1 0 1 2 0];
y = [16 9 17 12 22 13 8 15 19 11];

%% part 1-21

figure()
scatter(x,y, 90, 'fillcolor', 'k')
axis([-0.5 3.5 5 25])
title('Scatter Plot X vs Y', 'fontsize', 16)
xlabel('Number of Transfers', 'fontsize', 16)
ylabel('Number of Broken Ampules', 'fontsize', 16)


avg_x = mean(x)
avg_y = mean(y)

n = size(x,2);

Sxy = sum(x.*y) - (sum(x)*sum(y))/n
Sxx = sum(x.*x) - (sum(x)^2)/n

b1 = Sxy/Sxx

b0 = avg_y - b1*avg_x

%% part 2.6

% a - estimate B1 with 95% CI

y_hat = b1.*x + b0

e = y-y_hat

sum(e)

SSE = sum(e.^2)

prob = .975;

Zt = tinv(prob, 8)

Probt = tcdf(8.53, 8)