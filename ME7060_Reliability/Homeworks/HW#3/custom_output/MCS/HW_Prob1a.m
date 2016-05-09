%There are 3 random variables with following distributions and statistics
% - P is the concentrated load at center of beam, which is normally
%   distributed. These are the following statistics.
%   [50 kip, 10 kip]
%
% - E is the Young's Modulus of the beam, which is normally
%   distributed. These are the following statistics.
%   [29x10e6 lb/in^2, 1x10e5 lb/in^2]
%
% - w is the distributed load on the beam, which is log normally
%   distributed. These are the following statistics.
%   [1 kip/ft, 0.1 kip/ft]

clear all
clc

%Area Moment of Inertia
Inrt = 1.33*10^3; % in^4
%Convert length L into inches
L=30*12; % in

%Convert units to lb/in
mux_w = 83.33; % lb/in
sigx_w = 8.33; %lb/in

%Next convert log normal statistics mux and sigx into muy and sigy
sigy_w = log((sigx_w/mux_w)^2 + 1);
muy_w = log(mux_w) - 0.5*sigy_w^2;

%Define statistics for concentrated load
mux_P = 50000; %lbf
sigx_P = 10000; %lbf

%Define statistics for Young's Modulus
mux_E = 29*10^6; %lb/in^2
sigx_E = 1*10^5; %lb/in^2

%Following 2 lines define number of samples(N) and variables(NVar)
N = 1000000;
NVar = 3;

%Define random vector U(N, NVar) dimensions
U=rand(N,NVar);

%Next define distributions for variables P , E, and w
XP = mux_P + sigx_P*icdf('normal', U(:,1), 0, 1);

XE = mux_E + sigx_E*icdf('normal', U(:,2), 0, 1);

Xw = exp(muy_w + sigy_w*icdf('normal', U(:,3), 0, 1));

%Calculate deflection of midspan
def = (XP*L^3)./(48*XE*Inrt) + (5*Xw*L^4)./(385*XE*Inrt);

%Get statistics
mu_del = mean(def)
sig_del = std(def)

%Get failure probability
nf = sum(def>2.0)
pf=nf/N

%Store Statistics
fileID = fopen('C:\Users\WSUadm\Desktop\Reliability_HW3\custom_output\' ...
'MCS\MCSStats.txt', 'a');
fprintf(fileID, '%s %s %s \r\n', num2str(mu_del), num2str(sig_del), ...
    num2str(pf));
fclose(fileID);
%fclose(fid);