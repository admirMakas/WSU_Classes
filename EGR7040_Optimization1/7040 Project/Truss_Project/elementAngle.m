%CALCULATE ELEMENT ANGLE
function[el_Angle] = elementAngle(nodeDef, elementDef, numEl)
el_Angle = zeros(numEl, 1);
for i=1:numEl
    pointStart = elementDef(i, 1);
    pointEnd = elementDef (i, 2);
    x1 = nodeDef(pointStart, 1);
    y1 = nodeDef(pointStart, 2);
    x2 = nodeDef(pointEnd, 1);
    y2 = nodeDef(pointEnd, 2);
    el_Angle(i, 1) = atand((y2-y1)/(x2-x1)) ;
end