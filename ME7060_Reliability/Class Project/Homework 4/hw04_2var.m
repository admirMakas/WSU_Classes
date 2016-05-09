%% Homework 4 2variables                         Student:      Daniel Clark
clear all                                      % Instructor: Dr. Ha-Rok Bae
close all                                      % Class: ME 7060 Spring 2016
clc
format shorte
%% Problem 
% The failure event is defined as $delta_{max} > 2.0$. Anser the questions
% below.

% First Order Second Moment (FOSM, also known as MVFOSM method)

%% Displacement Equation
L = 30*12;                  % was ft now in
I = 1.33*10^3;              % in^4
W = 0.1/12;                 % was kip/ft now lbs/in
delta_max = @(P,E) (P.*L.^3)./ (48.*E.*I) + (5*W.*L.^4)./(385.*E.*I);
LSB_local = @(E) (2 - (5*W.*L.^4)./(385.*E.*I) ).* ( (48.*E.*I) / L.^3 );

mu_P = 50*10^3;             % was kip now lbs
mu_E = 29*10^6;             % lb/in^2
sigma_P = 10*10^3;          % was kip now lbs
sigma_E = 1*10^6;           % lb/in^2

S = [mu_E-3*sigma_E , mu_P-3*sigma_P; mu_E+3*sigma_E , mu_P+3*sigma_P];
gridNum = 20;
testpoints = gridsamp([min(S(:,1)) min(S(:,2)); max(S(:,1)) max(S(:,2))], gridNum);
Es = reshape(testpoints(:,1),gridNum,gridNum);
Ps = reshape(testpoints(:,2),gridNum,gridNum);
Responses = delta_max(testpoints(:,1),testpoints(:,2));
Rs = reshape(Responses,gridNum,gridNum);

figure,hold on
mesh(Es,Ps,Rs,40,'LineWidth',2)
break
Eplot_pts = linspace(mu_E-3*sigma_E,mu_E+3*sigma_E);
plot(Eplot_pts,LSB_local(Eplot_pts),'k','LineWidth',2)

ellipse(sigma_E,sigma_P,0,mu_E,mu_P,'k');
ellipse(2*sigma_E,2*sigma_P,0,mu_E,mu_P,'k');
ellipse(3*sigma_E,3*sigma_P,0,mu_E,mu_P,'k');

axis([min(S(:,1)) max(S(:,1)), min(S(:,2)) max(S(:,2))])
scatter(mu_E,mu_P,50,'o','MarkerFaceColor',[.7,.7,.7],'MarkerEdgeColor','k')
xlabel('E')
ylabel('P')
X1 = [mu_P; mu_E; mu_W];
[g, gDelta] = LSB(X1)
gtil = g + gDelta(1) * (X1(1) - mu_P) + gDelta(2) * (X1(2) - mu_E) ...
    + gDelta(3) * (X1(3) - mu_W)
sigma_g = sqrt( (gDelta(1) * X1(1))^2 + (gDelta(2) * X1(2))^2 ...
    + (gDelta(3) * X1(3))^2 )
beta = gtil/sigma_g

pf = 1 - cdf('normal',beta,0,1)


break


%% Constants
L = 30*12;                  % was ft now in
I = 1.33*10^3;              % in^4



%% Assumed Normal Variable
mu_W = 1/12;                % was kip/ft now lbs/in
sigma_W = 0.1/12;           % was kip/ft now lbs/in

%% Making arrays 
Number_of_runs = 100;
P_v = mu_P + sigma_P.*icdf('normal',rand(Number_of_runs,1),0,1);
E_v = mu_E + sigma_E.*icdf('normal',rand(Number_of_runs,1),0,1);
W_v = mu_W + sigma_W.*icdf('normal',rand(Number_of_runs,1),0,1);
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




%% log Normal
mu_w = 1*10^3/12;                   % was kip/ft now lbs/in
sigma_w = 0.1*10^3/12;              % was kip/ft now lbs/in
mu_log_w = log((mu_w^2)/sqrt(sigma_w^2+mu_w^2));
sigma_log_w = sqrt(log((sigma_w^2/(mu_w^2))+1));
PDF = @(x) exp(-1*((log(x) - mu_log_w).^2)/...
    (2*sigma_log_w^2)).*(1./(x.*sigma_log_w*sqrt(2*pi)));

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





