function [dN] = getdN(r, s, t) %Computes dN for stiffness matrix calc.
    
dN1r = -((s - 1)*(t - 1))/8;
dN1s = -((r - 1)*(t - 1))/8;
dN1t = -((r - 1)*(s - 1))/8;

dN2r = ((s - 1)*(t - 1))/8;
dN2s = ((r + 1)*(t - 1))/8;
dN2t = ((r + 1)*(s - 1))/8;

dN3r = -((s + 1)*(t - 1))/8;
dN3s = -((r + 1)*(t - 1))/8;
dN3t = -((r + 1)*(s + 1))/8;

dN4r = ((s + 1)*(t - 1))/8;
dN4s = ((r - 1)*(t - 1))/8;
dN4t = ((r - 1)*(s + 1))/8;

dN5r = ((s - 1)*(t + 1))/8;
dN5s = ((r - 1)*(t + 1))/8;
dN5t = ((r - 1)*(s - 1))/8;

dN6r = -((s - 1)*(t + 1))/8;
dN6s = -((r + 1)*(t + 1))/8;
dN6t = -((r + 1)*(s - 1))/8;

dN7r = ((s + 1)*(t + 1))/8;
dN7s = ((r + 1)*(t + 1))/8;
dN7t = ((r + 1)*(s + 1))/8;

dN8r = -((s + 1)*(t + 1))/8;
dN8s = -((r - 1)*(t + 1))/8;
dN8t = -((r - 1)*(s + 1))/8;

dN=[[dN1r dN1s dN1t]' [dN2r dN2s dN2t]' [dN3r dN3s dN3t]'...
    [dN4r dN4s dN4t]' [dN5r dN5s dN5t]' [dN6r dN6s dN6t]'...
    [dN7r dN7s dN7t]' [dN8r dN8s dN8t]'];

end