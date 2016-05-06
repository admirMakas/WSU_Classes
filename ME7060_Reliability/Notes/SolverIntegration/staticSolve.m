% clc;
% clear all;
% close all;
% format short g;
% ----------------------------------------------------------------------- %
% fprintf('#######################################################################\n')
% fprintf('############ A simple FEA solver for plane truss problems #############\n')
% fprintf('#### Developed by: Koorosh Gobal (gobal.2@wright.edu) - March 2016 ####\n')
% fprintf('#######################################################################\n')
% fprintf('\n');
function staticSolve
%% Read input file
nodeCoord = [];
elm = [];
area = [];
E = [];
dispBC = [];
forceBC = [];

fInput = fopen('input_file.inp');
line = fgetl(fInput);
while ~strcmp(line, 'END INPUT')
    line = fgetl(fInput);
    %% Read node coordinate
    if strcmp(line, 'BEGIN NODE_COORDINATE')
        while ~strcmp(line, 'END NODE_COORDINATE')
            line = fgetl(fInput);
            if strcmp(line, 'END NODE_COORDINATE')
                break
            end
            nodeCoord = [nodeCoord; str2num(line)];
        end
    end
    %% Read elements (node connectivity)
    if strcmp(line, 'BEGIN ELEMENTS')
        while ~strcmp(line, 'END ELEMENTS')
            line = fgetl(fInput);
            if strcmp(line, 'END ELEMENTS')
                break
            end
            elm = [elm; str2num(line)];
        end
    end
    %% Read elements (node connectivity)
    if strcmp(line, 'BEGIN ELEMENT_AREA')
        while ~strcmp(line, 'END ELEMENT_AREA')
            line = fgetl(fInput);
            if strcmp(line, 'END ELEMENT_AREA')
                break
            end
            area = [area; str2num(line)];
        end
    end
    %% Read modulus of elasticity
    if strcmp(line, 'BEGIN E')
        while ~strcmp(line, 'END E')
            line = fgetl(fInput);
            if strcmp(line, 'END E')
                break
            end
            E = [E; str2num(line)];
        end
    end
    %% Read displacement boundary
    if strcmp(line, 'BEGIN DISP_BC')
        while ~strcmp(line, 'END DISP_BC')
            line = fgetl(fInput);
            if strcmp(line, 'END DISP_BC')
                break
            end
            dispBC = [dispBC; str2num(line)];
        end
    end
    %% Read force boundary
    if strcmp(line, 'BEGIN F_BC')
        while ~strcmp(line, 'END F_BC')
            line = fgetl(fInput);
            if strcmp(line, 'END F_BC')
                break
            end
            forceBC = [forceBC; str2num(line)];
        end
    end
end
fclose(fInput);

nodeCoord = nodeCoord(:, 2:3);
elm = elm(:, 2:3);
area = area(:, 2);
E = E(:, 2);

%% Define degrees of freedom associated with each elemet
elmDOF = zeros(size(elm,1), 4);
for i = 1:size(elm,1)
    elmDOF(i, :) = [2 * elm(i, 1) - 1, 2 * elm(i, 1), 2 * elm(i, 2) - 1, 2 * elm(i, 2)];
end 
 %% Calculate element length
elmLength = zeros(size(elm,1), 1);
for i = 1:size(elm,1)
    elmLength(i) = ...
        sqrt((nodeCoord(elm(i, 2), 2) - nodeCoord(elm(i, 1), 2)).^2 + ...
             (nodeCoord(elm(i, 2), 1) - nodeCoord(elm(i, 1), 1)).^2);
end

%% Calculate element angle
elmAngle = zeros(size(elm,1), 1);
for i = 1:size(elm,1)
    elmAngle(i) = ...
        atan2(nodeCoord(elm(i, 2), 2) - nodeCoord(elm(i, 1), 2), ...
              nodeCoord(elm(i, 2), 1) - nodeCoord(elm(i, 1), 1));
end

%% Put stiffness matrix together
K = zeros(size(nodeCoord, 1) * 2);
for i = 1:size(elm, 1)
    theta = elmAngle(i);
    c = cos(theta); s = sin(theta);
%     rMat = [-c^2, c*s, -c^2, -c*s;
%             c*s, s^2, -c*s, -s^2;
%             -c^2, -c*s, c^2, c*s;
%             -c*s, -s^2, c*s, s^2];
    rMat = [c^2, c*s, -c^2, -c*s;
            c*s, s^2, -c*s, -s^2;
            -c^2, -c*s, c^2, c*s;
            -c*s, -s^2, c*s, s^2];
    k = E(i) * area(i) / elmLength(i) * rMat;
    K(elmDOF(i, :), elmDOF(i, :)) = K(elmDOF(i, :), elmDOF(i, :)) + k;
end

%% Assign boundary conditions
F = zeros(2 * size(nodeCoord, 1), 1);
for i = 1:size(forceBC, 1)
    F(2 * forceBC(i, 1) - 1) = forceBC(i, 2);
    F(2 * forceBC(i, 1)) = forceBC(i, 3);
end
remove = [];
for i = 1:size(dispBC, 1)
    if dispBC(i, 2)  && dispBC(i, 3)
        remove = [remove, 2 * dispBC(i, 1) - 1, 2 * dispBC(i, 1)];
    elseif dispBC(i, 2)
        remove = [remove, 2 * dispBC(i, 1) - 1];
    elseif dispBC(i, 3)
        remove = [remove, 2 * dispBC(i, 1)];
    end
end
K(remove, :) = [];
K(:, remove) = [];
F(remove, :) = [];
%% SOLVE
u = K \ F;

%% Write to output file
DOF = 1:2 * size(nodeCoord, 1);
U = zeros(length(DOF), 1);
U(setdiff(DOF, remove)) = u;
dlmwrite('U.out', U);
end