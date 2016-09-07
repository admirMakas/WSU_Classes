k=3;
n=3;
p=[0:0.001:1];
%p0=0.2155; .7840
%p1=0.1055; .8940
%p2=0.036; .9640
%p3=0.049; .9510

p = .9510;

bc = factorial(n)/(factorial(k)*factorial(n-k));

P = bc*p.^k.*(1-p).^(n-k);