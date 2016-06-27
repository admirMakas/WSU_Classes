%Desing variables are the following:
% x1 = outside diameter of the shaft, m
% x2 = ratio of inside/outside diameter (dimensionless)

function [] = solve(x)

clear all

% set options for the fmincon function

options = optimset('LargeScale', 'off', ...
    'TolCon', 1e-8, 'TolX', 1e-8, ...
    'Display', 'final-detailed', 'PlotFcns', @optimplotfval, ...
    'MaxFunEvals', 5000, 'MaxIter', 1000);

%Line of code below will display all the options for fmincon
%options = optimset('fmincon')

%Set bounds on the design variables
%Lb = [20; 0.60]; Ub = [500; 0.999];

%Set initial design
x0 = [10; 10];

% invoke fmincon function, four instances of "[]" indicate we have
% no linear constrains for this problem

[x, FunVal, ExitFlag, Output] = ...
    fminsearch(@ObjAndGrad, x0, options)

function [f, gf] = ObjAndGrad(x)

%rename design variables x
x1=x(1); x2=x(2);

%define const function
f = x1^2 + 10*x2^2 - 3*x1*x2;

%compute gradients of the objective function
%use nargout to determine if gf output is desired by user
if nargout > 1
    gf(1,1) = 2*x1 - 3*x2;
    gf(2,1) = 20*x2 - 3*x1;
end