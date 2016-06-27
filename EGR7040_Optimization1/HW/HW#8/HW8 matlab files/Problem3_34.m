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
Lb = [20; 0.60]; Ub = [500; 0.999];

%Set initial design
x0 = [80; 0.4];

% invoke fmincon function, four instances of "[]" indicate we have
% no linear constrains for this problem

[x, FunVal, ExitFlag, Output] = ...
    fmincon(@ObjAndGrad, x0, [], [], [], [], Lb, ...
    Ub, @ConstAndGrad, options)

function [f, gf] = ObjAndGrad(x)

%rename design variables x
x1=x(1); x2=x(2);

%define const function
f=(3.08269e-3)*(x1^2)*(1-x2^2);

%compute gradients of the objective function
%use nargout to determine if gf output is desired by user
if nargout > 1
    gf(1,1) = -(3554099593036481*x1*(x2^2 - 1))/576460752303423488;
    gf(2,1) = -(3554099593036481*x1^2*x2)/576460752303423488;
end

function [g, h, gg, gh] = ConstAndGrad(x)

% g returns inequality constraints
% h returns equality constraints
% gg returns gradients of the inequality constraints
% gh returns gradients of the equality constraints

% rename design variables
x1=x(1); x2=x(2);

% inequality constraints
g(1) = ((5.093e7)/(x1^3*(1-x2^4))) - 275;
g(2) = ((6.36619e5)/(x1^4*(1-x2^4))) - 3.49066e-2;
g(3) = 2.0e7 - 4.17246e4*x1^3*(1-x2)^(2.5);

%Problem has no equality constraints so we get
h=[];

%Compute constraint gradients
%use nargout to determine if gg and gh are requested by user
if nargout > 2
    gg(1,1) = 152790000/(x1^4*(x2^4 - 1));
    gg(2,1) = (203720000*x2^3)/(x1^3*(x2^4 - 1)^2);
    gg(1,2) = 2546476/(x1^5*(x2^4 - 1));
    gg(2,2) = (2546476*x2^3)/(x1^4*(x2^4 - 1)^2);
    gg(1,3) = -(17203756074113433*x1^2*(1 - x2)^(5/2))/137438953472;
    gg(2,3) = (28672926790189055*x1^3*(1 - x2)^(3/2))/274877906944;
    gh = [];
end
