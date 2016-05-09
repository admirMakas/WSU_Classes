close all
clear all
clc

%P E I
mu=[4 2*10^7 10^-4]
sig=[1 0.5*10^7 0.2*10^-4]
emu=mu;
esig=sig;
xk_iter=[];
uk_iter=[];
b_iter=[];
%k=1 %P E I
xk=mu

for i = 1:20
% Equiv norm variable for P
%extreme (Gumbel) dist parameters
alp=1.2825/sig(1)
d1=mu(1)-0.5772/alp
pk=xk(1);
fpk=alp*exp(-(pk-d1)*alp-exp(-(pk-d1)*alp))
Fpk=exp(-exp(-(pk-d1)*alp))
sig_equiv=pdf('norm',icdf('norm',Fpk,0,1),0,1)/fpk
mu_equiv=pk-icdf('norm',Fpk,0,1)*sig_equiv
emu(1)=mu_equiv;
esig(1)=sig_equiv;
uk=(xk-emu)./esig
xk_iter=[xk_iter;xk];
uk_iter=[uk_iter;uk];
%P E I
gxk=xk(2)*xk(3)-78.12*xk(1)
dgx=[-78.12, xk(3) , xk(2)]

B1=gxk-dgx(1)*esig(1)*uk(1)-dgx(2)*esig(2)*uk(2)-dgx(3)*esig(3)*uk(3)
B2=sqrt((dgx(1)*esig(1))^2+(dgx(2)*esig(2))^2+(dgx(3)*esig(3))^2)
bk=B1/B2
a1=-dgx(1)*esig(1)/B2
a2=-dgx(2)*esig(2)/B2
a3=-dgx(3)*esig(3)/B2
b_iter=[b_iter;bk]
xk_next=[emu(1)+bk*esig(1)*a1 emu(2)+bk*esig(2)*a2 emu(3)+bk*esig(3)*a3];
xk=xk_next

end


plot(b_iter,'*')
figure, scatter(xk_iter(:,1),xk_iter(:,2),'filled','MarkerFaceColor','k')


