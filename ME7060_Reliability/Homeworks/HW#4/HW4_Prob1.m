% Problem 1 code and result.

clear all
clc
Emu  = 29*10^6; % lb/in^2
Esig = 2*10^6; % lb/in^2

Pmu = 50000; % lb
Psig = 10000; % lb

Wmu = 83.33; % lb/in
Wsig = 8.33; % lb/in

L = 360; % in
I = 1330; % in^4

% First need to define the limit state function. We do not want beam
% displacement to be greater than 2 inches. Therefore following limit
% state function shall be used.
% -------------------------------------------------------------------------
% g(P, E, W) = 2 - [(PL^3)/(48EI) + (5/385)*((W*L^4)/(E*I))] = 0
% -------------------------------------------------------------------------
% If g(P, E, W)<0 then we have failure

% To determine failure probability 'Pf' need to get beta, which is
% beta = g_mu/g_sig

% Define g_mu
g_mu = 2-((Pmu*L^3)/(48*Emu*I) + (5/385)*((Wmu*L^4)/(Emu*I)));

% To calculate g_sig requires that the derivative of g(P, E, W) are taken
% with respect to all the random variable (P, E, W). Then need to take the
% magnutude of the derivatives in order to obtain g_sig.

% Get derivatives
syms P1 E1 W1

g = 2-((P1*L^3)/(48*E1*I) + (5/385)*((W1*L^4)/(E1*I)));

gP = diff(g, P1);
% display('deriv gP')
% eval(subs(gP, E1, Emu))
gE = diff(g, E1);
% display('deriv gE')
% eval(subs(gE, [P1, E1, W1], [Pmu, Emu, Wmu]))
gW = diff(g, W1);
% display('deriv gW')
% eval(subs(gW1, E1, Emu))

% Now get g_sig
g_sig = sqrt((eval(subs(gP, E1, Emu))*Psig)^2 + ...
    (eval(subs(gE, [P1, E1, W1], [Pmu, Emu, Wmu]))*Esig)^2 + ...
    (eval(subs(gW, E1, Emu))*Wsig)^2);

% Now get beta
beta = g_mu/g_sig

% Now get probability of failure
Pf = 1 - cdf('norm', beta, 0 ,1)