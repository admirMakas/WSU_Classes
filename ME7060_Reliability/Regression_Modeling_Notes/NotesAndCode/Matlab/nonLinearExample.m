%% Example nonlinear                             Student:      Daniel Clark
clear all                                      % Instructor: Dr. Ha-Rok Bae
close all                                      % Class: ME 7060 Spring 2016
clc
format shorte
%% Given
delta_max = @(P,L,E,I,w) (P.*L.^3)./ (48.*E.*I) + (5*w.*L.^4)./(385.*E.*I);
%constants
mu_E = 29*10^6;        % lb/in^2
mu_L = 30*12;          % was ft now in
mu_w = 1*10^3/12;      % was kip/ft now lbs/in

%% Normal
mu_I = 1.33*10^3;      % in^4
sigma_I = 9*10^1;      % in^4
mu_P = 50*10^3;        % was kip now lbs
sigma_P = 10*10^3;     % was kip now lbs

%% Lets look at the response surface
gridNum = 20;
S = [mu_P-3*sigma_P , mu_I-3*sigma_I; mu_P+3*sigma_P , mu_I+3*sigma_I];
testpoints = gridsamp([min(S(:,1)) min(S(:,2)); max(S(:,1)) max(S(:,2))], gridNum);
oneM = ones(1,gridNum*gridNum)';

M = [testpoints(:,1), oneM*mu_L, oneM*mu_E, testpoints(:,2),oneM*mu_w];

Y_true = delta_max(M(:,1),M(:,2),M(:,3),M(:,4),M(:,5));
Es = reshape(testpoints(:,1),gridNum,gridNum);
Is = reshape(testpoints(:,2),gridNum,gridNum);
Response = reshape(Y_true, size(Es));

Samples = [S(1,1),S(1,2); 
           S(2,1),S(1,2); 
           S(1,1),S(2,2);
           S(2,1),S(2,2)];
oneM = ones(1,length(Samples))';
M = [Samples(:,1), oneM*mu_L, oneM*mu_E, Samples(:,2), oneM*mu_w];
SampleResponse = delta_max(M(:,1),M(:,2),M(:,3),M(:,4),M(:,5));

%% Plot of true surface and 2^k points
figure, hold on
contour(Es,Is,Response,30);
scatter(Samples(:,1),Samples(:,2),100,'o','MarkerFaceColor','k','MarkerEdgeColor','k')
xlabel('P','FontSize',16);
ylabel('I','FontSize',16);

figure
mesh(Es,Is,Response)
xlabel('P','FontSize',16);
ylabel('I','FontSize',16);

%% Linear Regression
Tables = fitlm(Samples,SampleResponse);
betas = table2array(Tables.Coefficients(:,1));

Y_hat = betas(1) + testpoints(:,1)*betas(2) + testpoints(:,2)*betas(3);
Response_est = reshape(Y_hat, size(Es));

figure, hold on
surf(Es,Is,Response_est,'FaceColor','r','EdgeColor','k')
alpha(0.50)
mesh(Es,Is,Response)
scatter3(Samples(:,1),Samples(:,2),SampleResponse,100,'o','MarkerFaceColor','k','MarkerEdgeColor','k')
view(37.5*(-1),30)
xlabel('P','FontSize',16);
ylabel('I','FontSize',16);

%% Error plot
oneM = ones(1,1000)';
ppoints = linspace(mu_P-3*sigma_P, mu_P+3*sigma_P,1000)';
errorpoints = [ppoints,oneM*(mu_I+2*sigma_I)];

M = [errorpoints(:,1), oneM*mu_L, oneM*mu_E, errorpoints(:,2),oneM*mu_w];
Error_Response = delta_max(M(:,1),M(:,2),M(:,3),M(:,4),M(:,5));

Error_Y_hat = betas(1) + errorpoints(:,1)*betas(2) + errorpoints(:,2)*betas(3);

figure, hold on
plot(ppoints, Error_Response,'k','linewidth',2)
plot(ppoints, Error_Y_hat,'r','linewidth',2)
xlabel('P','FontSize',16);
ylabel('Max Deflection','FontSize',16);
legend('True Response','Linear Regression')

%% Nonlinear regression
Samples = [S(1,1),S(1,2); 
           S(2,1),S(1,2); 
           S(1,1),S(2,2);
           S(2,1),S(2,2);
           mu_P,mu_I];
oneM = ones(1,length(Samples))';
M = [Samples(:,1), oneM*mu_L, oneM*mu_E, Samples(:,2), oneM*mu_w];
SampleResponse = delta_max(M(:,1),M(:,2),M(:,3),M(:,4),M(:,5));

Tables = fitlm(Samples,SampleResponse,'interactions');
betas = table2array(Tables.Coefficients(:,1));
Y_hat = betas(1) + testpoints(:,1)*betas(2) + testpoints(:,2)*betas(3) + testpoints(:,1).*testpoints(:,2)*betas(4);
Response_est = reshape(Y_hat, size(Es));

figure, hold on
surf(Es,Is,Response_est,'FaceColor','r','EdgeColor','k')
alpha(0.50)
mesh(Es,Is,Response)
scatter3(Samples(:,1),Samples(:,2),SampleResponse,100,'o','MarkerFaceColor','k','MarkerEdgeColor','k')
view(37.5*(-1),30)
xlabel('P','FontSize',16);
ylabel('I','FontSize',16);

%% Error plot
Error_Y_hat = betas(1) + errorpoints(:,1)*betas(2) + errorpoints(:,2)*betas(3) + ...
    errorpoints(:,1).*errorpoints(:,2)*betas(4);

figure, hold on
plot(ppoints, Error_Response,'k','linewidth',2)
plot(ppoints, Error_Y_hat,'r','linewidth',2)
xlabel('P','FontSize',16);
ylabel('Max Deflection','FontSize',16);
legend('True Response','Linear Regression')



