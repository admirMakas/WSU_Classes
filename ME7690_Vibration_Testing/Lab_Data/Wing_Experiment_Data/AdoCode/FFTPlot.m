function [fout, Yout] = FFTPlot(x, dt);

% s = inputname(1)

% get data size
N = length(x);

% calculate time period
T=N*dt;

% calculate frequency step
df = 1/T;


X = fft(x, N)*dt;
%Y(1)=[];
if mod(N,2)==0
    X1=X(1:ceil(N/2),:);
else
    X1=X(1:ceil(N/2)-1,:);
end


%Plot Fourier Coefficients
% figure
% plot(X,'bo')
% title('Fourier Coefficients in the Complex Plane');
% xlabel('Real Axis');
% ylabel('Imaginary Axis');

%Define amplitude for frequency plot
power = real(X1.*conj(X1))/(N*dt)*2;
%power = abs(X1).^2;

%Define frequency range
Fend = df*(N/2-1);
f=(0:df:Fend)';

%Plot frequency vs power
figure
plot(f,power, 'linewidth', 2)
xlabel('Frequency')
ylabel('Amp')
title(['Frequency Plot '])
axis([5,60,0,0.001]);
grid on

%return
if nargout>0
    fout = f;
    Xout = X1;
end
end