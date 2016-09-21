clc

data1 = [216 162 153 216 225 216 303 225 243 189];
data2 = [225 171 198 189 189 135 162 135 117 162];

avg1 = mean(data1)
avg2 = mean(data2)

stdev1 = std(data1)^2
stdev2 = std(data2)^2

stdevH = (9*stdev1 + 9*stdev2)/(18)

prob = .975;

Zt = tinv(prob, 9)

Probt = tcdf(1.237, 7)