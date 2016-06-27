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
%Lb = [20; 0.60]; Ub = [500; 0.999];

%Set initial design
x0 = [10; 10];

% invoke fmincon function, four instances of "[]" indicate we have
% no linear constrains for this problem

[x, FunVal, ExitFlag, Output] = ...
    fmincon(@ObjAndGrad, x0, [], [], [], [], [], ...
    [], @ConstAndGrad, options)

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

function [g, h, gg, gh] = ConstAndGrad(x)

% g returns inequality constraints
% h returns equality constraints
% gg returns gradients of the inequality constraints
% gh returns gradients of the equality constraints

% rename design variables
x1=x(1); x2=x(2);

% inequality constraints
%g(1) = -x1 - x2 - 5;
%g(2) = -x1 - 5;
%g(3) = x2 - 5;
%g(4) = -2*x1 - x2 + 4;

g=[];

%Problem has no equality constraints so we get
h=[];

%Compute constraint gradients
%use nargout to determine if gg and gh are requested by user
if nargout > 2
    %gg(1,1) = -1;
    %gg(2,1) = -1;
    %gg(1,2) = -1;
    %gg(2,2) = 0;
    %gg(1,3) = 1;
    %gg(2,3) = 0;
    %gg(1,4) = -2;
    %gg(2,4) = -1;
    gg=[];
    gh = [];
end
