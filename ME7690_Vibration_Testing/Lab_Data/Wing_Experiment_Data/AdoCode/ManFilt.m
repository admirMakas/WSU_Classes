function y = ManFilt(x, dt)

%get length of x
N=length(x);

%get period T
T=N*dt;

%get time vector
TimeRange=(0:dt:T)';
TimeRange=TimeRange(1:end-1,1);

%get frequency step
df=1/T;

%put time data into frequency domain
X=fft(x);

%define variable "del", which determines how many Fourier coefficients
%will be set to zero. For now the value will be kept constant 7.
del = 20;

%set first "del" values to zero in X
X(1:del,1) = 0;

yout=real(ifft(X));
y=yout;
