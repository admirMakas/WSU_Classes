%%
% This is the Newmark Method implementation
%

Mred=[1 0; 0 2];
Kred=[6 -2; -2 8];

b=0.25;
y=0.50;
dt = 0.24216;

Dn = zeros(2,1);
dDn = zeros(2,1);
ddDn = zeros(2,1);

Dn_Hist = [];
dDn_Hist = [];
ddDn_Hist = [];

Dn_Hist=[Dn_Hist, Dn];
dDn_Hist=[dDn_Hist, dDn];
ddDn_Hist=[ddDn_Hist, ddDn];

F0 = zeros(2,1);

t=0;
t_arr=[];
t_arr=[t_arr, t];
count=0;
for i=1:11
%while t<20 && isnan(Dn(1,1))~=1
   count = count+1;
   t=t+dt;
   t_arr=[t_arr, t];
   
   %if t<=0.01
       F0(2,1) = 10;
   %else
   %    F0(149,1) = 0;
   %end
   
   Dn1 = ( ( (1/(b*dt^2))*Mred + Kred )\eye(2) )* ...
       (F0 + Mred*( (1/(b*dt^2))*Dn + (1/(b*dt))*dDn + (1/(2*b)-1)*ddDn ));
   
   ddDn1 = (1/(b*dt^2))*( Dn1 - Dn - dt*dDn ) - (1/(2*b) - 1)*ddDn;
   
   dDn1 = (y/(b*dt))*(Dn1 - Dn) - (y/b - 1)*dDn - dt*(y/(2*b) - 1)*ddDn;
   
   Dn_Hist=[Dn_Hist, Dn1];
   dDn_Hist=[dDn_Hist, dDn1];
   ddDn_Hist=[ddDn_Hist, ddDn1];
   
   Dn = Dn1;
   dDn = dDn1;
   ddDn = ddDn1;
   
end

figure
plot(t_arr, Dn_Hist(1,:))
figure
plot(t_arr, dDn_Hist(1,:))
figure
plot(t_arr, ddDn_Hist(1,:))