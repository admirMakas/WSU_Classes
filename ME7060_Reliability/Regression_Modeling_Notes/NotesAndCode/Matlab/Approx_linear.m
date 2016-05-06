function [ g_l ] = Approx_linear(g0, x1, x0, dg_dxi)
g_l = g0;
for i = 1:length(x0)
g_l = g_l + ( x1(i) - x0(i) )*dg_dxi(i);
end
end