ans = zeros(1,350);

for R = 1:350
    syms t
    f=(4.04488e6)/(R^3*t + (R*t^3)/4) + (5257.21*(2*R+t))/(R^3*t + (R*t^3)/4) + (1.6492e10*(2*R + t))/(R^3*t + (R*t^3)/4)^2 - 1;
    s=solve(f,t);
    ans(R)=s(1);
end
ans';