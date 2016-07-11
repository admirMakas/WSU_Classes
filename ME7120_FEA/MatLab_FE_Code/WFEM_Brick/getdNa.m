function [dNa] = getdNa(r, s, t) %Computes dN for stiffness matrix calc.
    
% Shape function for incompatible modes:
dN9r = -2*r;
dN9s = 0;
dN9t = 0;

dN10r = 0;
dN10s = -2*s;
dN10t = 0;

dN11r = 0;
dN11s = 0;
dN11t = -2*t;

dNa=[[dN9r dN9s dN9t]'...
    [dN10r dN10s dN10t]'...
    [dN11r dN11s dN11t]'];

end