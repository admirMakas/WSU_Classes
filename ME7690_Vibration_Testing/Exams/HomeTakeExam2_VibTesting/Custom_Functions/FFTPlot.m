function [fout, Yout] = FFTPlot(x, dt);

% get data size
N = length(x);

% calculate time period
T=N*dt;

TimeRange=(0:dt:T)';
TimeRange=TimeRange(1:end-1,1);

% calculate frequency step
df = 1/T;

X = fft(x, N);
%Y(1)=[];
if mod(N,2)==0
    X1=X(1:ceil(N/2),:);
else
    X1=X(1:ceil(N/2)-1,:);
end


%Plot Fourier Coefficients
% figure
% plot(X(1:(N/2)),'bo')
title('Fourier Coefficients in the Complex Plane');
xlabel('Real Axis');
ylabel('Imaginary Axis');

%Define amplitude for frequency plot
power = real(X1.*conj(X1))/(N*dt)*2;
%pow = power(2:4096,1);
%power = abs(X1).^2;

%Define frequency range
Fend = df*(N/2-1);
f=(0:df:Fend)';

%Plot frequency vs power
% figure
% plot(f(1:N/2,1), power(1:N/2,1), 'r', 'linewidth', 2)
xlabel('Frequency')
title('Frequency Plot')

%Define inverse FFT for modified frequency data to see impact on time
%domain
del=10;
Xs=X;
Xs(1:del,1) = 0;
if mod(N,2)==0
    X1s=Xs(1:ceil(N/2),:);
else
    X1s=Xs(1:ceil(N/2)-1,:);
end
powers=real(X1s.*conj(X1s))/(N*dt)*2;
plot(f, powers, 'r', 'linewidth', 2)
hold on
plot(f, power)

y=ifft(Xs);
figure
plot(TimeRange, x)
hold on
plot(TimeRange, real(y))

%return
if nargout>0
    fout = f;
    Xout = X1;
end
end