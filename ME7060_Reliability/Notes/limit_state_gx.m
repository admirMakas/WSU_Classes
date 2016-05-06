function [gxk,dgxk,d2gxk]=limit_state_gx(x)

gxk=x(1)^4+2*x(2)^4-20;

dgxk=[4*x(1)^3 2*4*x(2)^3]';
d2gxk=[12*x(1)^2 2*12*x(2)^2]';



% % 
% % P=x(1);
% % E=x(2);
% % I=x(3);
% % gxk=E*I-78.12*P;




