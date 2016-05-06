function [ outputText ] = figureTag(figNumber)

st1 = '00';
st2 = num2str(figNumber);
st1(end-length(st2)+1:end) = [];
outputText = [st1,st2];

end