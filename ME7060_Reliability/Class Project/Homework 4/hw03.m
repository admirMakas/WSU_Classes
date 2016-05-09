%% Homework 3                                    Student:      Daniel Clark
clear all                                      % Instructor: Dr. Ha-Rok Bae
close all                                      % Class: ME 7060 Spring 2016
clc
format shorte
% Given
delta_max = @(P,L,E,I,w) (P.*L.^3)./ (48.*E.*I) + (5*w.*L.^4)./(385.*E.*I);
L = 30*12;                  % was ft now in
I = 1.33*10^3;              % in^4

%% Normal
mu_P = 50*10^3;             % was kip now lbs
mu_E = 29*10^6;             % lb/in^2
sigma_P = 10*10^3;          % was kip now lbs
sigma_E = 1*10^5;           % lb/in^2

%% log Normal
mu_w = 1*10^3/12;                   % was kip/ft now lbs/in
sigma_w = 0.1*10^3/12;              % was kip/ft now lbs/in
mu_log_w = log((mu_w^2)/sqrt(sigma_w^2+mu_w^2));
sigma_log_w = sqrt(log((sigma_w^2/(mu_w^2))+1));
PDF = @(x) exp(-1*((log(x) - mu_log_w).^2)/...
    (2*sigma_log_w^2)).*(1./(x.*sigma_log_w*sqrt(2*pi)));

%% Making arrays 
Number_of_runs = 100;
P_v = mu_P + sigma_P.*icdf('normal',rand(Number_of_runs,1),0,1);
E_v = mu_E + sigma_E.*icdf('normal',rand(Number_of_runs,1),0,1);
W_v = exp(mu_log_w + sigma_log_w.*icdf('normal',rand(Number_of_runs,1),0,1));
L_v = ones(Number_of_runs,1)*L;
I_v = ones(Number_of_runs,1)*I;

%% Over-plot data on PDFs
figure, hold on
xP = linspace(30*10^3,75*10^3,1000);
P = normpdf(xP,mu_P,sigma_P);
plot(xP,P)
scatter(P_v, normpdf(P_v,mu_P,sigma_P),100,'.')

figure, hold on
xE = linspace(28*10^6,30*10^6,1000);
E = normpdf(xE,mu_E,sigma_E);
plot(xE,E)
scatter(E_v, normpdf(E_v,mu_E,sigma_E),100,'.')

figure, hold on
xW = linspace(60,120,1000);
W = PDF(xW);
plot(xW,W)
scatter(W_v, PDF(W_v),100,'.')

%% Stastics 100 runs
responseVector = delta_max(P_v,L_v,E_v,I_v,W_v);
MEAN_response_100 = mean(responseVector)
STDEV_response_100 = std(responseVector)

%% Stastics 10000 runs
Number_of_runs = 10000;
P_v = mu_P + sigma_P.*icdf('normal',rand(Number_of_runs,1),0,1);
E_v = mu_E + sigma_E.*icdf('normal',rand(Number_of_runs,1),0,1);
W_v = exp(mu_log_w + sigma_log_w.*icdf('normal',rand(Number_of_runs,1),0,1));
L_v = ones(Number_of_runs,1)*L;
I_v = ones(Number_of_runs,1)*I;

responseVector = delta_max(P_v,L_v,E_v,I_v,W_v);
MEAN_response_10000 = mean(responseVector)
STDEV_response_10000 = std(responseVector)

%% Stastics 10^6 runs
Number_of_runs = 10^6;
P_v = mu_P + sigma_P.*icdf('normal',rand(Number_of_runs,1),0,1);
E_v = mu_E + sigma_E.*icdf('normal',rand(Number_of_runs,1),0,1);
W_v = exp(mu_log_w + sigma_log_w.*icdf('normal',rand(Number_of_runs,1),0,1));
L_v = ones(Number_of_runs,1)*L;
I_v = ones(Number_of_runs,1)*I;

responseVector = delta_max(P_v,L_v,E_v,I_v,W_v);
MEAN_response_10_6 = mean(responseVector)
STDEV_response_10_6 = std(responseVector)

%% LHS design 

Number_of_runs = 10000;
LHSArray = lhsdesign(Number_of_runs,3,'iterations',10);

P_v = mu_P + sigma_P.*icdf('normal',LHSArray(:,1),0,1);
E_v = mu_E + sigma_E.*icdf('normal',LHSArray(:,2),0,1);
W_v = exp(mu_log_w + sigma_log_w.*icdf('normal',LHSArray(:,3),0,1));
L_v = ones(Number_of_runs,1)*L;
I_v = ones(Number_of_runs,1)*I;

responseVector = delta_max(P_v,L_v,E_v,I_v,W_v);
MEAN_response_LHS_10000 = mean(responseVector)
STDEV_response_LHS_10000 = std(responseVector)

figure, hold on
xP = linspace(0,100,1000);
P = normpdf(xP,mu_P,sigma_P);
plot(xP,P)
scatter(P_v, normpdf(P_v,mu_P,sigma_P),100,'.')

figure, hold on
xE = linspace(28*10^6,30*10^6,1000);
E = normpdf(xE,mu_E,sigma_E);
plot(xE,E)
scatter(E_v, normpdf(E_v,mu_E,sigma_E),100,'.')

figure, hold on
xW = linspace(0,2,1000);
W = PDF(xW);
plot(xW,W)
scatter(W_v, PDF(W_v),100,'.')






break
%% convergence plot
runs = linspace(100,10^6,100)';
meanResponse = runs;
stdevResponse = runs;

tic
parfor i = 1:length(runs)
    Number_of_runs = runs(i);
    P_v = mu_P + sigma_P.*icdf('normal',rand(Number_of_runs,1),0,1);
    E_v = mu_E + sigma_E.*icdf('normal',rand(Number_of_runs,1),0,1);
    W_v = exp(mu_log_w + sigma_log_w.*icdf('normal',rand(Number_of_runs,1),0,1));
    L_v = ones(Number_of_runs,1)*L;
    I_v = ones(Number_of_runs,1)*I;
    responseVector = delta_max(P_v,L_v,E_v,I_v,W_v);
    meanResponse(i) = mean(responseVector);
    stdevResponse(i) = std(responseVector);
end
toc

figure
plot(runs,meanResponse)

figure
plot(runs,stdevResponse)





