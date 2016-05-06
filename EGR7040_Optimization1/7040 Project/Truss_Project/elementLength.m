%CALCULATE ELEMENT LENGTH
function [el_Length] = elementLength(nodeDef, elementDef, numEl)
el_Length = zeros(numEl, 1);
for i =1:numEl
    pointStart = elementDef(i, 1);
    pointEnd = elementDef(i, 2);
    x1 = nodeDef(pointStart, 1);
    y1 = nodeDef(pointStart, 2);
    x2 = nodeDef(pointEnd, 1);
    y2 = nodeDef(pointEnd, 2);
    el_Length(i, 1) = sqrt(power(x2-x1, 2) + power(y2-y1, 2));
end