%Desing variables are the following:
% x1 = inner radius R, cm
% x2 = wall thickness t, cm

function [] = solve

clear all

% set options for the fmincon function

options = optimset('LargeScale', 'off', 'GradObj', 'on', ...
    'GradConstr', 'on', 'TolCon', 1e-8, 'TolX', 1e-8, ...
    'Display', 'final-detailed', 'PlotFcns', @optimplotfval, ...
    'MaxFunEvals', 5000, 'MaxIter', 1000, 'Algorithm', 'interior-point');

%Line of code below will display all the options for fmincon
%options = optimset('fmincon')

%Set bounds on the design variables
Lb = [35; 1]; Ub = [200; 20];

%Set initial design
x0 = [250; 15];

% Linear inequality constraints
A=[-2, -1; 1, -0.5; -1, 0.5; 0, 1; 0, -1];
b=[-70; 250; -35; 20; -1];
%A=[1, -0.5; -1, 0.5];
%b=[250; -35];


% invoke fmincon function, four instances of "[]" indicate we have
% no linear constrains for this problem

[x, FunVal, ExitFlag, Output] = ...
    fmincon(@ObjAndGrad, x0, [], [], [], [], [], ...
    [], @ConstAndGrad, options)

function [f, gf] = ObjAndGrad(x)

%rename design variables x
R=x(1); t=x(2);

%define const function
f = 153.71*R*t;

%compute gradients of the objective function
%use nargout to determine if gf output is desired by user
if nargout > 1
    gf(1,1) = 153.71*t;
    gf(2,1) = 153.71*R;
end

function [g, h, gg, gh] = ConstAndGrad(x)

% g returns inequality constraints
% h returns equality constraints
% gg returns gradients of the inequality constraints
% gh returns gradients of the equality constraints

% rename design variables
R=x(1); t=x(2);

% inequality constraints
g(1) = (4.04488e6)/(R^3*t + (R*t^3)/4) + (5257.21*(2*R + t))/(R^3*t + (R*t^3)/4) + (1.6492e10*(2*R + t))/(R^3*t + (R*t^3)/4)^2 - 1;
g(2) = 70 - (2*R + t);
g(3) = (2*R)/t - 91;
g(4) = 4.479e7/(R^3*t + (R*t^3)/4) - 20;
g(5) = R - 0.5*t - 250;
g(6) = 35 - (R + 0.5*t);
g(7) = t - 40;
g(8) = 1 - t;

%Problem has no equality constraints so we get
h=[];

%Compute constraint gradients
%use nargout to determine if gg and gh are requested by user
if nargout > 2
    gg(1,1) = 5780363524660265/(549755813888*(R^3*t + (R*t^3)/4)) - (4044880*(3*R^2*t + t^3/4))/(R^3*t + (R*t^3)/4)^2 + 32984000000/(R^3*t + (R*t^3)/4)^2 - (2*(3*R^2*t + t^3/4)*(32984000000*R + 16492000000*t))/(R^3*t + (R*t^3)/4)^3 - ((3*R^2*t + t^3/4)*((5780363524660265*R)/549755813888 + (5780363524660265*t)/1099511627776))/(R^3*t + (R*t^3)/4)^2;
    gg(2,1) = 780363524660265/(1099511627776*(R^3*t + (R*t^3)/4)) - (4044880*(R^3 + (3*R*t^2)/4))/(R^3*t + (R*t^3)/4)^2 + 16492000000/(R^3*t + (R*t^3)/4)^2 - (2*(32984000000*R + 16492000000*t)*(R^3 + (3*R*t^2)/4))/(R^3*t + (R*t^3)/4)^3 - (((5780363524660265*R)/549755813888 + (5780363524660265*t)/1099511627776)*(R^3 + (3*R*t^2)/4))/(R^3*t + (R*t^3)/4)^2;
    gg(1,2) = -2;
    gg(2,2) = -1;
    gg(1,3) = 2/t;
    gg(2,3) = -(2*R)/t^2;
    gg(1,4) = -(44790000*(3*R^2*t + t^3/4))/(R^3*t + (R*t^3)/4)^2;
    gg(2,4) = -(44790000*(R^3 + (3*R*t^2)/4))/(R^3*t + (R*t^3)/4)^2;
    gg(1,5) = 1;
    gg(2,5) = -0.5;
    gg(1,6) = -1;
    gg(2,6) = -0.5;
    gg(1,7) = 0;
    gg(2,7) = 1;
    gg(1,8) = 0;
    gg(2,8) = -1;
    gh = [];
end
