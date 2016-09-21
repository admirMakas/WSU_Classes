clc

Joggers = [20 19 22 21 53 22 20 31 28];

FPE = [79 27 54 28 31 42 37 24 23 22 32 23 25 41 98];

data = [Joggers FPE];

index = [1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2];

boxplot(data, index)

y1 = mean(Joggers)

y2 = mean(FPE)

var1 = std(Joggers)^2

var2 = std(FPE)^2

t_obs = (y1-y2)/sqrt(var1/9 + var2/15)

df = (var1/9 + var2/15)^2/((116.9/9)^2/8 + (494.5/15)^2/14)

prob = .975;

Zt = tinv(prob, 21)

Prob_t = tcdf(-1.8946, df)