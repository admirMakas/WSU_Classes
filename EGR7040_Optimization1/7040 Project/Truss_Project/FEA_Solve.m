
function [output] = FEA_Solve(area, output_request)
% DEFINIG NODES AND ELEMENTS
nodeDef = [1.5 2.598; 4.5 2.598; 7.5 2.598; 3 0; 6 0; 0 0; 9 0]; %[5 by 2] matrix
elementDef = [1 2; 2 3; 6 4; 4 5; 5 7; 6 1; 1 4; 4 2; 2 5; 5 3; 3 7]; %[6 by 2] matrix 

%area = 0.02.*ones(11,1); %initial area definition [10 by 1]

% MECHANICAL PROPERTIES OF ELEMENTS
youngMS = 200e9;
DOF = 2;

% DEFINE NODAL BOUNDARY CONDITIONS
% DISPALCEMENT( deltaU) = [NODE_NUMBER X_DISPLACEMENT Y_DISPLACEMENT]
% FORCE = [NODE_NUMBER F_X F_Y]
DispBC = [6 0 0; 7 0 0]; %displacement BCs [2 by 3]
force = [4 0 -5000000; 5 0 -5000000]; %force BCs [2 by 3]
numNode = size(nodeDef, 1); %scalar = 6
numEl = size(elementDef, 1); %scalar = 11
youngM = youngMS.*ones(numEl, 1); %Young's Modulus defined for all 10 elements [10 by 1]

% CALCULATING LENGTH OF EACH ELEMENT AND ALSO ITS ANGLE TO HORIZION
elLength = elementLength(nodeDef, elementDef, numEl);
elAngle = elementAngle(nodeDef, elementDef, numEl);

% CALCULATE GLOBAL STIFFNESS MATRIX
gsMat = GSM(elementDef, numNode, numEl, youngM, area, elLength, elAngle, DOF, DispBC);
% generates [8 by 8] matrix

% PUT THE FORCE MATRIX IN GLOBAL FORM
F = [];
for i = 1:size(force, 1) %goes from 1 to 2 since there are 2 loads applied
    F = [F; force(i, 2:3)']; %makes column force vector for containing applied load
end
tMat = zeros(size(gsMat, 1), size(F, 1)); %creates [8 by 4] zero matrix
for i = 1:size(force, 1) %For loop ensures that load values are in correct position withing the force matrix
    tMat(force(i ,1)*2 - 1, 2*i - 1) = 1;
    tMat(force(i ,1)*2, 2*i) = 1;
end
F = tMat*F; %creates [8 by 1] vector 

% DISPLACEMENT SOLVER
nodeDisp = inv(gsMat)*F;
U = nodeDisp;
nodeDisp = [nodeDisp; zeros(4, 1)];

% STRESS SOLVER
elStress = feaStress(elementDef, numEl, nodeDisp, elLength, youngM, area, elAngle);
if(strcmp(output_request, 'stress'))
    output = elStress;
elseif(strcmp(output_request, 'displacement'))
    output = nodeDisp;
end