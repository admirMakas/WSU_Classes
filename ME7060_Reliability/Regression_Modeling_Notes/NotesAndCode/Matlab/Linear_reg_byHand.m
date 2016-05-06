%% Example 0                                     Student:      Daniel Clark
clear all                                      % Instructor: Dr. Ha-Rok Bae
close all                                      % Class: ME 7060 Spring 2016
clc
format shorte
%% 
S = [1;2;3;4]

Y = [4;7;10;13]

%% linear 
F_lin = [ones(length(S),1),S]
beta = (F_lin'*F_lin)\(F_lin'*Y)
Tables_lm = fitlm(S,Y);

figure, hold on
scatter(S,Y,100,'o','MarkerFaceColor','k','MarkerEdgeColor','k')
x = linspace(0,4,100)';
y_hat = beta(1) + x*beta(2);
plot(x,y_hat,'r','linewidth',2)

%% nonlinear 
F_lin = [ones(length(S),1),S,S.^2]
beta = (F_lin'*F_lin)\(F_lin'*Y)
Tables_nlm = fitlm(S,Y,'purequadratic');

y_hat = beta(1) + x*beta(2) + x.^2*beta(3);
plot(x,y_hat,'b','linewidth',2)