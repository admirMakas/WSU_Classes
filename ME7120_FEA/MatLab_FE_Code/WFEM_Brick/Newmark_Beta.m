%% Used to get K and M matrices

% num_nodes = size(KModel,1)/6;
% keep_dof = [(0:(num_nodes-1))*6 + 1];
% keep_dof = sort([keep_dof, keep_dof + 2, keep_dof + 4]);
% keep_dof= [keep_dof(3:151), keep_dof(153)];
% Kred = KModel(keep_dof, keep_dof);
% Mred = Mmodel(keep_dof, keep_dof);

%%
% This is the Newmark Method implementation
%
b=0.25;
y=0.50;
dt = 0.0001;

Dn = zeros(150,1);
dDn = zeros(150,1);
ddDn = zeros(150,1);

Dn_Hist = [];
dDn_Hist = [];
ddDn_Hist = [];

Dn_Hist=[Dn_Hist, Dn];
dDn_Hist=[dDn_Hist, dDn];
ddDn_Hist=[ddDn_Hist, ddDn];

F0 = zeros(150,1);

t=0;
t_arr=[];
t_arr=[t_arr, t];
count=0;
while t<0.125 && isnan(Dn(1,1))~=1
   count = count+1;
   t=t+dt;
   t_arr=[t_arr, t];
   
   if t<=0.01
       F0(149,1) = 100000;
   else
       F0(149,1) = 0;
   end
   
   Dn1 = ( ( (1/(b*dt^2))*Mred + Kred )\eye(150) )* ...
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

% figure
% plot(t_arr, -Dn_Hist(121,:))
% figure
% plot(t_arr, -dDn_Hist(121,:))
% figure
% plot(t_arr, -ddDn_Hist(121,:))