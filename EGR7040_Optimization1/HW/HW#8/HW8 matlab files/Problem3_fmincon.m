%Desing variables are the following:
% x1 = outside diameter of the shaft, m
% x2 = ratio of inside/outside diameter (dimensionless)

function [] = solve(x)

clear all

% set options for the fmincon function

options = optimset('LargeScale', 'off', 'GradObj', 'on', ...
    'GradConstr', 'on', 'TolCon', 1e-8, 'TolX', 1e-8, ...
    'Display', 'final-detailed', 'PlotFcns', @optimplotfval, ...
    'MaxFunEvals', 5000, 'MaxIter', 1000, 'Algorithm', 'interior-point');

%Line of code below will display all the options for fmincon
%options = optimset('fmincon')

%Set bounds on the design variables
Lb = [0; 0]; Ub = [70; 50];

%Set initial design
x0 = [1; 1];

% invoke fmincon function, four instances of "[]" indicate we have
% no linear constrains for this problem

[x, FunVal, ExitFlag, Output] = ...
    fmincon(@ObjAndGrad, x0, [], [], [], [], Lb, ...
    Ub, @ConstAndGrad, options)

function [f, gf] = ObjAndGrad(x)

%rename design variables x
x1=x(1); x2=x(2);

%define const function
f=20*x1 + 64*x2;

%compute gradients of the objective function
%use nargout to determine if gf output is desired by user
if nargout > 1
    gf(1,1) = 20;
    gf(2,1) = 64;
end

function [g, h, gg, gh] = ConstAndGrad(x)

% g returns inequality constraints
% h returns equality constraints
% gg returns gradients of the inequality constraints
% gh returns gradients of the equality constraints

% rename design variables
x1=x(1); x2=x(2);

% inequality constraints
g(1) = -25*x1 - 70*x2 + 2100;

%Problem has no equality constraints so we get
h=[];

%Compute constraint gradients
%use nargout to determine if gg and gh are requested by user
if nargout > 2
    gg(1,1) = -25;
    gg(2,1) = -70;
    gh = [];
end
