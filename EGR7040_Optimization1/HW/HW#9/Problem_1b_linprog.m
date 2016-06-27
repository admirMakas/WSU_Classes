%Desing variables are the following:
% x1 = outside diameter of the shaft, m
% x2 = ratio of inside/outside diameter (dimensionless)

function [] = solve(x)

clear all

% set options for the fmincon function

options = optimset('LargeScale', 'off', ...
    'TolCon', 1e-8, 'TolX', 1e-8, ...
    'Display', 'final-detailed', 'PlotFcns', @optimplotfval, ...
    'MaxFunEvals', 5000, 'MaxIter', 1000, 'Algorithm', 'interior-point');

%Line of code below will display all the options for fmincon
%options = optimset('fmincon')

%Set bounds on the design variables
Lb = [0, 0];
Ub = [];
%Set initial design
x0 = [];

f=[8 -3 15];

A=[-5 1.8 3.6; -3 -6 -8.2; -1.5 4 -7.5; 0 1 -5]

b=[-2 -5 4.5 -1.5];

% invoke fmincon function, four instances of "[]" indicate we have
% no linear constrains for this problem
%[x, FunVal, ExitFlag, Output]
[x, FunVal, ExitFlag, Output] = ...
    linprog(f, A, b, [], [], Lb, Ub, x0, options)