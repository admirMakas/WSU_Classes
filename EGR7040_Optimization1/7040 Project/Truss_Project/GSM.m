%ASSEMBLE GLOBAL STIFFNESS MATRIX
function globalStiffnessMatrix = GSM(elementDef, numNode, numEl, youngM, area, elLength, elAngle, DOF, DispBC)
globalStiffnessMatrix = zeros(numNode * DOF, numNode * DOF); %matrix of zeros [12 by 12]
elStifness = zeros(numEl, 1); %element stiffness matrix initiated [10 by 1]
for i=1:numEl
    elStifness(i) = youngM(i) * area (i) / elLength(i);
end
for i=1:numEl %individual element global stiffness matrices calculated and added together to construct total global stiffness matrix [14 by 14]
    phi = elAngle(i);
    l2gMat = zeros(numNode * DOF, numNode * DOF);
    nodeStart = elementDef(i, 1);
    nodeEnd = elementDef(i, 2);
    
    rMat = [cosd(phi) sind(phi) 0 0 ; ...
    -sind(phi) cosd(phi) 0 0; ...
    0 0 cosd(phi) sind(phi); ...
    0 0 -sind(phi) cosd(phi)]; %rotation matrix defined here
    elsMat = rMat' * (elStifness(i) * [1 0 -1 0; 0 0 0 0; -1 0 1 0; 0 0 0 0]) * rMat; %rotated element stiffness matrix, total of 10 of them [4 by 4]

    %4 lines below add a "1" in the appropriate location within 12gMat,
    %this is used to generate a [14 by 14] matrix correspondig to a
    %specific truss element
    l2gMat(2 * nodeStart - 1, 2 * nodeStart - 1) = 1;
    l2gMat(2 * nodeStart, 2 * nodeStart) = 1;
    l2gMat(2 * nodeEnd - 1, 2 * nodeEnd - 1) = 1;
    l2gMat(2 * nodeEnd, 2 * nodeEnd) = 1;
    
    %reconstructs 12gMat by keeping only non-zero columns
    l2gMat(:, ~any(l2gMat, 1)) = []; %columns
    
    elsMat = l2gMat * elsMat * l2gMat'; %Finally generate [12 by 12] stiffness matrix for subject truss element
    globalStiffnessMatrix = globalStiffnessMatrix + elsMat; %add up all the individual global truss element stiffness matrices to create global stiffness matrix
end

tMat = eye(size(globalStiffnessMatrix, 1)); %create [124 by 14] identity matrix
for i =1:size(DispBC, 1)
    tMat(DispBC(i, 1)*2-1, DispBC(i, 1)*2-1) = 0;
    tMat(DispBC(i, 1)*2, DispBC(i, 1) * 2) = 0;
end
tMat(:, ~any(tMat, 1)) = [];
globalStiffnessMatrix = tMat' * globalStiffnessMatrix * tMat;