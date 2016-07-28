% We don't have a lot of time to get fancy. Let's just use Guyan to
% remove DOFs 7-42, and the Lagrange 385-420 DOFs. Alternatively,
% we could try ditching the first 6 as well.
clear all
clc

load mk.mat
M11=M([43:384],[43:384]);
K11=K([43:384],[43:384]);
M12=M([1:6,7:42,385:420],[43:384]);
K12=K([1:6,7:42,385:420],[43:384]);
M21=M12';
K21=K12';
M22=M([1:6,7:42,385:420],[1:6,7:42,385:420]);
K22=K([1:6,7:42,385:420],[1:6,7:42,385:420]);
T=[eye(size(K11)); -K22\K21'];
Mred=T'*[M11 M21;M12 M22]*T;
Kred=T'*[K11 K21;K12 K22]*T;
X1=rand(size(Mred,1),20)-.5;
%Shift of 1.
shift=1;
Kredp=Kred+shift*Mred;
KinvM=real(Kredp\Mred);
%not required
%%%%%%%%%%%%%%
error=1;
i=0;
%%%%%%%%%%%%%%%%
% I did this with a for loop and "watched" during my "exam conditions"
while error>.00001
%%%%%%%%%%%%
i=i+1;
%%%%%%%%%%%%%%
X2=real(KinvM*X1);
Ksmall=X2'*Kredp*X2;
Msmall=X2'*Mred*X2;
disp('Msmall eig')
min(eig(Msmall))
[v,d]=eig(Msmall\Ksmall);
X1=X2*v;
freqs=sort(sqrt(real(diag(d))-shift)/2/pi)
%not required
if i>1
i
size(freqs)
size(oldfreqs)
error=abs((freqs(11)-oldfreqs(11))/freqs(11))
end
oldfreqs=freqs;
size(freqs)
size(oldfreqs)
%%%%%%%%%%%%
X1=X1/norm(X1);
%If you don't normalize somehow, it breaks down due to
%ill-conditioning pretty quickly (but not before you have two
%places in the eigenvalues.
end
%X1'*Kred*X1
max(abs(X1'*Kred*X1-diag(X1'*Kred*X1)))
max(abs(X1'*Mred*X1-diag(X1'*Mred*X1)))