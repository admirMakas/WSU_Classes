
function [] = solve(x)

clear;
clc;
%format short g;

area = 0.02.*ones(11,1); %initial area definition [11 by 1]

Lb = 0.0001 * ones (11,1); %lower bound on x-sec area
Ub = 0.070 * ones (11,1); %upper bound on x-sec area

options=optimset('Algorithm', 'interior-point', 'MaxIter', 1E5, 'MaxFunEvals', 1E5, 'TolCon', 1E5, ...
    'PlotFcns', @optimplotfval);

x0 = area; %initial truss element x-sec areas [11 by 1]

%fmincon call for optimization
[x, fval] = fmincon(@CostFunc, x0, [], [], [], [], Lb, Ub, @ConsFun, options)

%additional outputs
[area, x] %vector output of original area and optimized areas [11 by 2]
[FEA_Solve(area, 'stress'), FEA_Solve(x, 'stress')] %vector output of original truss stress and optimized stress [11 by 2]
[FEA_Solve(area, 'displacement'), FEA_Solve(x, 'displacement')] %vector output of original truss displacement and optimized displacement [11 by 2]
[CostFunc(area), CostFunc(x)] %vector output of starting and optimized costFun values
[ConsFun(area), ConsFun(x)] %vector output of stress constraints for each truss 1-11=tension 12-22=compression [22 by 2]

%COST FUNCTION DEFINITION
function f = CostFunc(x)
density = 7850; %material density definition
L = [3.0; 3.0; 3.0; 3.0; 3.0; 3.0; 3.0; 3.0; ...
    3.0; 3.0; 3.0]; %truss element lenghts
f = density * (x(1)*L(1) + x(2)*L(2) + x(3)*L(3) + x(4)*L(4) + x(5)*L(5) + ...
x(6)*L(6) + x(7)*L(7) + x(8)*L(8) + x(9)*L(9) + x(10)*L(10) + x(11)*L(11)); %cost function defined in terms of design variables

%CONSTRAINT FUNCTION DEFINITION, NO GRADIENTS SPECIFIED
function [g, h] = ConsFun(area)

stress = FEA_Solve(area, 'stress'); %element stresses [11 by 1]
sigma_max = 550e6;
g = [stress./sigma_max-1; -stress./sigma_max - 1]; %inequality constraint definitions
h = [];
