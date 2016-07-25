% close all
% clear all
% clc
% 
% syms r s
% 
% N1 = (1-r-s)*(1-2*r-2*s);
% N2 = r*(2*r-1);
% N3 = s*(2*s-1);
% N4 = 4*r*(1-r-s);
% N5 = 4*r*s;
% N6 = 4*s*(1-r-s);
% 
% x = [0, 1, 0, 0.5, 0.5, 0]';
% y = [0, 0, 1, 0, 0.5, 0.75]';
% 
% N = [N1;N2;N3;N4;N5;N6]
% [sum(diff(N,r).*x), sum(diff(N,r).*y);
% sum(diff(N,s).*x), sum(diff(N,s).*y)]

syms e

dN = [0.5*(-1+2*e) -2*e 0.5*(1+2*e)]
x=[0 0.75 1.0]

%e=1

%eval(dN)*x'

(dN)*x';

%.5*(-1+2*e) * (-2*e)


(-2/(e-1))*(-2*e*(e - 1/2))*(0.5-e/2)