%% Single Point Approx                           Student:      Daniel Clark
clear all                                      % Instructor: Dr. Ha-Rok Bae
close all                                      % Class: ME 7060 Spring 2016
clc
format shorte
%% Given information
E = 29*10^6;        % lb/in^2
L = 30*12;          % was ft now in
w = 1*10^3/12;      % was kip/ft now lbs/in

delta_max = @(P,I) (P.*L.^3)./ (48.*E.*I) + (5*w.*L.^4)./(385.*E.*I);
dmax_dP = @(P,I) (L.^3)./ (48.*E.*I);
dmax_dI = @(P,I) (-P.*L.^3)./ (48.*E.*I.^2) - (5*w.*L.^4)./(385.*E.*I.^2);

%% Normal
mu_I = 1.33*10^3;      % in^4
sigma_I = 9*10^1;      % in^4
mu_P = 50*10^3;        % was kip now lbs
sigma_P = 10*10^3;     % was kip now lbs

x2 = [mu_P, mu_I];
g1 = delta_max(mu_P, mu_I);

dg_dxi = [dmax_dP(mu_P, mu_I),dmax_dI(mu_P, mu_I)];

%% Compare the approximation Methods
gridNum = 20;
S = [mu_P-3*sigma_P , mu_I-3*sigma_I; mu_P+3*sigma_P , mu_I+3*sigma_I];
testpoints = gridsamp([min(S(:,1)) min(S(:,2)); max(S(:,1)) max(S(:,2))], gridNum);


number = gridNum*gridNum;
matrix = zeros(number,3);
direction = [1,1,1];
StepSize = linspace(-5,5,number);

for i = 1:number
    xApprox = testpoints(i,:);
    act = delta_max(testpoints(i,1),testpoints(i,2));
    L = Approx_linear(g1, xApprox, x2, dg_dxi);
    R = Approx_reciprocal(g1, xApprox, x2, dg_dxi);
    
    matrix(i,1) = L;
    matrix(i,2) = R;
    matrix(i,3) = act;

end

%% Plots
Ps = reshape(testpoints(:,1),gridNum,gridNum);
Is = reshape(testpoints(:,2),gridNum,gridNum);
Linear_approx = reshape(matrix(:,1), size(Ps));
Reciprocal_approx = reshape(matrix(:,2), size(Ps));
true = reshape(matrix(:,3), size(Ps));

figure, hold on
surf(Ps,Is,Linear_approx,'FaceColor','r','EdgeColor','k')
mesh(Ps,Is,true)
view(37.5*(-1),30)
xlabel('P','FontSize',16);
ylabel('I','FontSize',16);

figure, hold on
surf(Ps,Is,Reciprocal_approx,'FaceColor','r','EdgeColor','k')
mesh(Ps,Is,true)
view(37.5*(-1),30)
xlabel('P','FontSize',16);
ylabel('I','FontSize',16);



