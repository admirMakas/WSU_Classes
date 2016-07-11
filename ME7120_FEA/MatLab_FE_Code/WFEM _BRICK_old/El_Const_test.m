clear all
clc

num_gauss=2;
[int_p, int_w] = gauss(num_gauss);
intPts=zeros(num_gauss^3,3);
intWts=zeros(num_gauss^3,3);
index=0;

for i=1:num_gauss
    for j=1:num_gauss
        for k=1:num_gauss
            index=index+1;
            intPts(index,:) = [int_p(i) int_p(j) int_p(k)];
            intWts(index,:) = [int_w(i) int_w(j) int_w(k)];
        end
    end
end

E = 2.03e11; % Pa
v = 0.29;
rho = 8050; % kg

%define global element coordinates (in this case same as local coordinates)

%Natural Coordinates used for visual check
r=[-1; 1; 1; -1; -1; 1; 1; -1];
s=[-1; -1; 1; 1; -1; -1; 1; 1];
t=[-1; -1; -1; -1; 1; 1; 1; 1];

X=1*[-1; 1; 1; -1; -1; 1; 1; -1];
Y=1*[-1; -1; 1; 1; -1; -1; 1; 1];
Z=1*[-1; -1; -1; -1; 1; 1; 1; 1];

scatter3(r, s, t, 'filled')
hold on
scatter3(X, Y, Z, 'filled')
hold on
scatter3(intPts(:,1), intPts(:,2), intPts(:,3), 'filled')

global_nodes = [[X(1) Y(1) Z(1)]; [X(2) Y(2) Z(2)]; [X(3) Y(3) Z(3)]; [X(4) Y(4) Z(4)];...
    [X(5) Y(5) Z(5)]; [X(6) Y(6) Z(6)]; [X(7) Y(7) Z(7)]; [X(8) Y(8) Z(8)]];

B=zeros(6,24);
Ke=zeros(24,24);
Me=zeros(24,24);
id=[1 4 7 10 13 16 19 22];

Em=getE(E,v);

% Loop to construct BRICK stiffness matrix
for p=1:num_gauss^3
    r = intPts(p,1);
    s = intPts(p,2);
    t = intPts(p,3);
    
    Ne=getN(r,s,t);
    dN=getdN(r,s,t);
    
    J=dN*global_nodes;
    JDet = det(J);
    
    for q=1:8
        dN_i = dN(:,q);
        
        Bi=[dot(J(1,:),dN_i) 0 0;
            0 dot(J(2,:),dN_i) 0;
            0 0 dot(J(3,:),dN_i);
            dot(J(2,:),dN_i) dot(J(1,:),dN_i) 0;
            0 dot(J(3,:),dN_i) dot(J(2,:),dN_i);
            dot(J(3,:),dN_i) 0 dot(J(1,:),dN_i)];
        B(1:end, id(q):id(q)+2) = Bi(1:end, 1:end);
    end
    
    Ki = prod(intWts(p,1:end))*JDet*(B'*Em*B);
    Ke=Ke+Ki;
    
    Mi = rho*prod(intWts(p,1:end))*JDet*(Ne'*Ne);
    Me=Me+Mi;
end
