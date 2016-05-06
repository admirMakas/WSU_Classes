function [ g_r ] = Approx_reciprocal(g0, x1, x0, dg_dxi)
g_r = g0;
for i = 1:length(x0)
g_r = g_r + ( x1(i) - x0(i) ) * ( x0(i) / x1(i) ) * dg_dxi(i);
end
end