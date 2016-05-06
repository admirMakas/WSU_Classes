function xp = derivFUN_vdp(t,x)
xp = zeros(2,1);
xp(1) = x(2);
xp(2) = x(2)*(1-x(1)^2)-x(1) + cos(t)+2*cos(5*t);