function [guk,dguk,d2guk]=limit_state_gu(u)

mu=[10 10]';
sig=[5 5]';

x=u(:).*sig+mu;

[gxk,dgxk,d2gxk]=limit_state_gx(x);

guk=gxk;
dguk=dgxk.*sig;
d2guk=d2gxk.*sig.^2;

% % Transformation
% %
% P=u(1);
% E=u(2);
% I=u(3);
% 
% 
% mu=[4 2*10^7 10^-4]';
% sig=[1 0.5*10^7 0.2*10^-4]';
% 
% x=u(:).*sig+mu;

% 
% % when u(1)=P~Type I extreme value Gumble(mx, sx)
% % mup=delta+0.5772/alpha
% % sigp=1.2825/alpha
% %
% alpha=1.2825
% delta=3.5499
% PHI_u=cdf('norm',u(1),0,1)
% x(1)=(alpha*delta-log(-log(PHI_u)))/alpha; % inverse extreme type-I pp.105~106
% 

% X-space function evaluation
% [gxk]=limit_state_gx(x);


