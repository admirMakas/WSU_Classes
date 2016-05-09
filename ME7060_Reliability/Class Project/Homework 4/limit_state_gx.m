function [gxk, dgxk]=limit_state_gx(x)

% gxk=x(1)^3+x(2)^3-18;



%
P=x(1);
E=x(2);
I=x(3);
gxk=E*I-78.12*P;




