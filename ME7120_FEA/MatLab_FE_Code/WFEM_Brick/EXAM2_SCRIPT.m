clear all
clc
%%
%Code below will generate stiffness matrix for 3 noded bar element
%user can supply global coordinates and code takes care of rest

num_gauss=3;
[int_p, int_w] = gauss(num_gauss)
%   intPts=zeros(num_gauss^3,3);
%   intWts=zeros(num_gauss^3,3);

syms e

N = [0.5*(-e+e^2) 1-e^2 0.5*(e+e^2)];

dN = diff(N,e);

x=[0 0.5 1];

J=dN*x';

K=zeros(3,3);
M=zeros(3,3);
K2 = zeros(3,3);
M2 = zeros(3,3);
for p=1:num_gauss
    
    e=int_p(p);
    
    Ni = eval(N);
    Bi=eval(dN);
    Ji=eval(J);
    
    Ki = int_w(p)*(Bi'*Bi)*(1/Ji);
    Mi = int_w(p)*(Ni'*Ni)*(Ji);
    
%     Ki2 = int_w(p)*(Bi'*Bi);
%     Mi2 = int_w(p)*(Ni'*Ni);
    
    K=K+Ki;
    M=M+Mi;
    
%     K2 = K2 + Ki2;
%     M2 = M2 + Mi2;
    
end

K
M

% K2
% M2

%%
% %Code below will plot Jacobian for Q4 element
% 
% a=[0 1 2 4];
% n=-1:0.1:1;
% 
% for l=1:length(a)
%     x=[1 0 -1 0;
%         0 a(l) 0 -1];
%     
%     J = @(e,n) 0.25*[-(1-n) 1-n 1+n -(1+n); -(1-e) -(1+e) 1+e 1-e]*x';
%     
%     Jdet=[];
%     for m = 1:length(n)
%         
%         Jdet_i=det(J(n(m),-n(m)));
%         Jdet=[Jdet Jdet_i];
%         
%     end
%     
%     plot(n,Jdet)
%     hold on
%     
% end

%%
%Code below will generate stiffness matrix for 2 noded beam element
%user can supply global coordinates and code takes care of rest

% num_gauss=2;
% [int_p, int_w] = gauss(num_gauss);
% %   intPts=zeros(num_gauss^3,3);
% %   intWts=zeros(num_gauss^3,3);
% 
% syms e
% 
% % N = [0.25*e^3 - 3/4*e + 0.5 0.25*(e^3 - e^2 - e + 1) ...
% %     0.25*(e^3 + 3*e + 2) 0.25*(e^3 + e^2 - e - 1)];
% 
% N = [0.25*(1-e)^2*(2+e) (1/8)*(1-e)^2*(1+e) ...
%     0.25*(1+e)^2*(2-e) (1/8)*(1+e)^2*(1-e)];
% 
% dN = diff(N,e);
% dN2 = diff(N,e,2)
% 
% % x=[0 1];
% 
% % J=[dN2(1) dN2(2)]*x';
% 
% J=0.5;
% 
% K=zeros(4,4);
% M=zeros(4,4);
% K2 = zeros(4,4);
% M2 = zeros(4,4);
% for p=1:num_gauss
%     
%     e=int_p(p);
%     
%     Ni = eval(N);
%     Bi=eval(dN2);
%     %Ji=eval(J);
%     
%     Ki = int_w(p)*(Bi'*Bi)*(1/Ji);
%     Mi = int_w(p)*(Ni'*Ni)*(Ji);
%     
% %     Ki2 = int_w(p)*(Bi'*Bi);
% %     Mi2 = int_w(p)*(Ni'*Ni);
%     
%     K=K+Ki;
%     M=M+Mi;
%     
% %     K2 = K2 + Ki2;
% %     M2 = M2 + Mi2;
%     
% end
% 
% K
% M
% 
% % K2
% % M2
% 