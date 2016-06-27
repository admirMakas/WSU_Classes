function [] = solve(x)

clear all

% set options for the fmincon function

options = optimset('LargeScale', 'off', 'GradObj', 'on', ...
    'GradConstr', 'on', 'TolCon', 1e-8, 'TolX', 1e-8, ...
    'Display', 'iter', 'PlotFcns', @optimplotfval);

%Line of code below will display all the options for fmincon
%options = optimset('fmincon')

%Set bounds on the design variables
Lb = [13; 0]; Ub = [100; 100];

%Set initial design
x0 = [20.1; 5.84];

% invoke fmincon function, four instances of "[]" indicate we have
% no linear constrains for this problem

[x, FunVal, ExitFlag, Output] = ...
    fmincon(@ObjAndGrad7_3, x0, [], [], [], [], Lb, ...
    Ub, @ConstAndGrad7_3, options)

function [f, gf] = ObjAndGrad7_3(x)

%rename design variables x
x1=x(1); x2=x(2);

%define const function
f=(x1 - 10)^3 + (x2 - 20)^3;

%compute gradients of the objective function
%use nargout to determine if gf output is desired by user
if nargout > 1
    gf(1,1) = 3*(x1 - 10)^2;
    gf(2,1) = 3*(x2 - 20)^2;
end

function [g, h, gg, gh] = ConstAndGrad7_3(x)

% g returns inequality constraints
% h returns equality constraints
% gg returns gradients of the inequality constraints
% gh returns gradients of the equality constraints

% rename design variables
x1=x(1); x2=x(2);

% inequality constraints
g(1) = 100 - (x1-5)^2 - (x2 - 5)^2;
g(2) = -82.81 + (x1 - 6)^2 + (x2 - 5)^2;

%Problem has no equality constraints so we get
h=[];

%Compute constraint gradients
%use nargout to determine if gg and gh are requested by user
if nargout > 2
    gg(1,1) = -2*(x1 - 5);
    gg(2,1) = -2*(x2 - 5);
    gg(1,2) = 2*(x1 - 6);
    gg(2,2) = 2*(x2 - 5);
    gh = [];
end
