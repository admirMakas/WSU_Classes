% Linear tapered beam
numberOfNodes = 33;

length = 25; % zero to this value
tbegin = 0.02; % initial width
tend = 0.01; % final width

width = @(x) (tend-tbegin)/length * x + tbegin;
Ixx = @(x)  1/12*width(x)^4;
Iyy = @(x)  1/12*width(x)^4;
J = @(x) 0.95*(Ixx(x)+Iyy(x));

width(length)
no = linspace(0,length,numberOfNodes)';
numberOfElements = numberOfNodes - 1;
% Info = zeros(numberOfElements,8);
InfoMat = cell(numberOfElements,9);
for i = 1:numberOfElements % A1     A3   J1 J3 Ixx1 Ixx3 Iyy1 Iyy3
    Info = [ width(no(i))^2,  width(no(i+1))^2, J(no(i)), ...
        J(no(i+1)), Ixx(no(i)), Ixx(no(i+1)), Iyy(no(i)), Iyy(no(i+1))];
    
    InfoMat{i,1} = 'steel';
    
    for j = 1:8
        InfoMat{i,1+j} = Info(j);
    end
end

nodes = [[1:numberOfNodes]',zeros(numberOfNodes,1), no, zeros(numberOfNodes,1)]


n1n2pnP = zeros(numberOfElements,4);

for i = 1:numberOfElements
    n1n2pnP(i,1) = i;
    n1n2pnP(i,2) = i+1;
    n1n2pnP(i,3) = 1;
    n1n2pnP(i,4) = i; % n1n2pnP(i,4) = i;
end
    n1n2pnP
    
    
    


