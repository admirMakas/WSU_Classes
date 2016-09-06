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
clear all
clc

tic

load Proj3_MK.mat

b = 1/4;
y = 1/2;
dt = 0.0001;
time_lim = 0.125;

Isp = eye(size(Kred,1));

num_col = int64(ceil(time_lim/dt));

Dn = zeros(150,1);
dDn = zeros(150,1);
ddDn = zeros(150,1);

Dn_Hist = zeros(150,num_col);
dDn_Hist = zeros(150,num_col);
ddDn_Hist = zeros(150,num_col);

Dn_Hist(:,1) = Dn;
dDn_Hist(:,1) = dDn;
ddDn_Hist(:,1) = ddDn;

F0 = zeros(150,1);

t=0;
t_arr = zeros(1,num_col);
t_arr(1,1) = t;
count = 1;

for i=1:num_col
    if isnan(Dn(1,1))~=1
        
        count = count+1;
        t = t+dt;
        t_arr(1,count) = t;
        
        if t<=0.01
            F0(149,1) = 100000;
        else
            F0(149,1) = 0;
        end
        
        Dn1 = ( ( (1/(b*dt^2))*Mred + Kred )\Isp )* ...
            (F0 + Mred*( (1/(b*dt^2))*Dn + (1/(b*dt))*dDn + (1/(2*b)-1)*ddDn ));
        
        ddDn1 = (1/(b*dt^2))*( Dn1 - Dn - dt*dDn ) - (1/(2*b) - 1)*ddDn;
        
        dDn1 = (y/(b*dt))*(Dn1 - Dn) - (y/b - 1)*dDn - dt*(y/(2*b) - 1)*ddDn;
        
        Dn_Hist(:,count) = Dn1;
        dDn_Hist(:,count) = dDn1;
        ddDn_Hist(:,count) = ddDn1;
        
        Dn = Dn1;
        dDn = dDn1;
        ddDn = ddDn1;
        
    else
        
        break
        
    end
end

figure
plot(t_arr, -Dn_Hist(121,:))
figure
plot(t_arr, -dDn_Hist(121,:))
figure
plot(t_arr, -ddDn_Hist(121,:))

toc