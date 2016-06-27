clear all
ans = zeros(1,350);

for R = 1:350
    syms t
    f=(4.479e7)/(R^3*t + (R*t^3)/4) - 20;
    s=solve(f,t);
    ans(R)=s(1);
end
ans';