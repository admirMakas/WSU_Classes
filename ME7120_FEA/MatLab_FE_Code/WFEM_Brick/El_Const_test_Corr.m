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

Bd=zeros(6,24);
Ba=zeros(6,9);
Me=zeros(24,24);

Kbb = zeros(24,24);
Kba = zeros(24,9);
Kab = zeros(9,24);
Kaa = zeros(9,9);

Em=getE(E,v);

%==========================================================================
%Need dN evaluated at location r=s=t=0 for the incompatible modes
dN0=getdN(0,0,0);

%Get Jacobian evaluated at the center of the brick element used for
%Gauss integration of incompatible modes
J0=dN0*global_nodes;
Jinv0=J0\eye(3);
%==========================================================================

% Loop to construct BRICK stiffness matrix
for p=1:num_gauss^3
    r = intPts(p,1);
    s = intPts(p,2);
    t = intPts(p,3);
    
    Ne=getN(r,s,t);
    dN=getdN(r,s,t);
    
    %Derivatives of incompatible shape functions
    dNa=getdNa(r,s,t);
    
    J=dN*global_nodes;
    Jinv=J\eye(3);
    JDet = det(J);
    
    for q=1:11
        if q<=8
            dN_i = dN(:,q);
            
            Bi=[Jinv(1,:)*dN_i 0 0;
                0 Jinv(2,:)*dN_i 0;
                0 0 Jinv(3,:)*dN_i;
                Jinv(2,:)*dN_i Jinv(1,:)*dN_i 0;
                0 Jinv(3,:)*dN_i Jinv(2,:)*dN_i;
                Jinv(3,:)*dN_i 0 Jinv(1,:)*dN_i];
            
            Bd(1:end, 1+(q-1)*3:1+(q-1)*3+2) = Bi(1:end, 1:end);
        else
            dN_i = dNa(:,q-8);
            
            Bi=[Jinv0(1,:)*dN_i 0 0;
                0 Jinv0(2,:)*dN_i 0;
                0 0 Jinv0(3,:)*dN_i;
                Jinv0(2,:)*dN_i Jinv0(1,:)*dN_i 0;
                0 Jinv0(3,:)*dN_i Jinv0(2,:)*dN_i;
                Jinv0(3,:)*dN_i 0 Jinv0(1,:)*dN_i];
            
            Ba(1:end, 1+(q-1-8)*3:1+(q-1-8)*3+2) = Bi(1:end, 1:end);
        end
    end
    
    Kbbi = prod(intWts(p,1:end))*JDet*(Bd'*Em*Bd);
    Kbai = prod(intWts(p,1:end))*JDet*(Bd'*Em*Ba);
    Kabi = prod(intWts(p,1:end))*JDet*(Ba'*Em*Bd);
    Kaai = prod(intWts(p,1:end))*JDet*(Ba'*Em*Ba);
    
    Kbb = Kbb+Kbbi;
    Kba = Kba+Kbai;
    Kab = Kab+Kabi;
    Kaa = Kaa+Kaai;
    
    Mi = rho*prod(intWts(p,1:end))*JDet*(Ne'*Ne);
    Me=Me+Mi;
end

Ke = Kbb - Kba*(Kaa\Kab);
Kg=zeros(72,72);

nodes = [1 2 3 4 5 6 7 8;
    5 6 7 8 9 10 11 12];
%nodes = [5 6 7 8 9 10 11 12];
%indices = zeros(1,24);
for i=1:2
    indices = zeros(1,24);
    for w = 1:8
        
        node_i = nodes(i,w);
        
        indices(3*w-2:3*w) = 1+(node_i-1)*6:3+(node_i-1)*6;
        
    end
    Kg(indices,indices) = Kg(indices,indices) + Ke;
end