clear all
clc

load mk.mat
M11=M([1:6,43:384],[1:6, 43:384]);
K11=K([1:6,43:384],[1:6, 43:384]);
M12=M([7:42,385:420],[1:6,43:384]);
K12=K([7:42,385:420],[1:6,43:384]);
M21=M12';
K21=K12';
M22=M([7:42,385:420],[7:42,385:420]);
K22=K([7:42,385:420],[7:42,385:420]);

%solve for the eigenvalues using eig in order to compare to the sub-space
%iteration 

[~,e] = eig(full(M11)\full(K11));
FreqCheck = sort(sqrt(real(diag(e))))/2/pi;

Q=[eye(size(K11)); -K22\K21'];
Mred=Q'*[M11 M21;M12 M22]*Q;
Kred=Q'*[K11 K21;K12 K22]*Q;
X1=rand(size(Kred,1),20)-.5;
%Shift of 1.
shift=1;
Kshift=Kred+shift*Mred;
KinvM=real(Kshift\Mred);
%not required
%%%%%%%%%%%%%%
error=1;
i=0;
%%%%%%%%%%%%%%%%
while error>.001
%%%%%%%%%%%%
i=i+1;
%%%%%%%%%%%%%%
X2=real(KinvM*X1);
Ksmall=X2'*Kshift*X2;
Msmall=X2'*Mred*X2;
disp('Msmall eig')
min(eig(Msmall))
[v,d]=eig(Msmall\Ksmall);
X1=X2*v;
freqs=sort(sqrt(real(diag(d))-shift)/2/pi)
%not required
% if i>1
% i
% size(freqs)
% size(oldfreqs)
% error=abs((freqs(11)-oldfreqs(11))/freqs(11))
% end
oldfreqs=freqs;
size(freqs)
size(oldfreqs)
%%%%%%%%%%%%
%Normalize
X1=X1/norm(X1);
end
%X1'*Kred*X1
max(abs(X1'*Kred*X1-diag(X1'*Kred*X1)))
max(abs(X1'*Mred*X1-diag(X1'*Mred*X1)))