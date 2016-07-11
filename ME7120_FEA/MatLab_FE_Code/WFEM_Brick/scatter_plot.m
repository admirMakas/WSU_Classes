% clear all
% clc
% 
% r=[-1; 1; 1; -1; -1; 1; 1; -1];
% s=[-1; -1; 1; 1; -1; -1; 1; 1];
% t=[-1; -1; -1; -1; 1; 1; 1; 1];
% 
% % scatter3(dr2, ds2, dt2, 'filled')
% num_gauss=3
% [int_p, int_w] = gauss(num_gauss)
% intPts=zeros(num_gauss^3,3);
% intWts=zeros(num_gauss^3,3);
% index=0;
% 
% for i=1:num_gauss
%     for j=1:num_gauss
%         for k=1:num_gauss
%             index=index+1
%             intPts(index,:) = [int_p(i) int_p(j) int_p(k)];
%             intWts(index,:) = [int_w(i) int_w(j) int_w(k)];
%         end
%     end
% end
% 
% scatter3(r, s, t, 'filled')
% hold on
% scatter3(intPts(:,1), intPts(:,2), intPts(:,3),'filled')

clc

syms r s t real
% 
% N1 = 1/8*((1-r)*(1-s)*(1-t));
% N2 = 1/8*((1+r)*(1-s)*(1-t));
% N3 = 1/8*((1+r)*(1+s)*(1-t));
% N4 = 1/8*((1-r)*(1+s)*(1-t));
% 
% N5 = 1/8*((1-r)*(1-s)*(1+t));
% N6 = 1/8*((1+r)*(1-s)*(1+t));
% N7 = 1/8*((1+r)*(1+s)*(1+t));
% N8 = 1/8*((1-r)*(1+s)*(1+t));
% 
% diff(N8, r)
% diff(N8, s)
% diff(N8, t)


dN1r = -((s - 1)*(t - 1))/8;
dN1s = -((r - 1)*(t - 1))/8;
dN1t = -((r - 1)*(s - 1))/8;

dN2r = ((s - 1)*(t - 1))/8;
dN2s = ((r + 1)*(t - 1))/8;
dN2t = ((r + 1)*(s - 1))/8;

dN3r = -((s + 1)*(t - 1))/8;
dN3s = -((r + 1)*(t - 1))/8;
dN3t = -((r + 1)*(s + 1))/8;

dN4r = ((s + 1)*(t - 1))/8;
dN4s = ((r - 1)*(t - 1))/8;
dN4t = ((r - 1)*(s + 1))/8;

dN5r = ((s - 1)*(t + 1))/8;
dN5s = ((r - 1)*(t + 1))/8;
dN5t = ((r - 1)*(s - 1))/8;

dN6r = -((s - 1)*(t + 1))/8;
dN6s = -((r + 1)*(t + 1))/8;
dN6t = -((r + 1)*(s - 1))/8;

dN7r = ((s + 1)*(t + 1))/8;
dN7s = ((r + 1)*(t + 1))/8;
dN7t = ((r + 1)*(s + 1))/8;

dN8r = -((s + 1)*(t + 1))/8;
dN8s = -((r - 1)*(t + 1))/8;
dN8t = -((r - 1)*(s + 1))/8;

dN=[[dN1r dN1s dN1t]' [dN2r dN2s dN2t]' [dN3r dN3s dN3t]'...
    [dN4r dN4s dN4t]' [dN5r dN5s dN5t]' [dN6r dN6s dN6t]'...
    [dN7r dN7s dN7t]' [dN8r dN8s dN8t]']