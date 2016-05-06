% x = input
% y = output
function [Gxx] = AutoSD(x, y, dt)

% get data size
N = length(x);

% determine # of columns
Col = size(x);

% calculate time period
T = N*dt;

% calculate frequency step
df = 1/T;

% put data into frequency domain
% for both x and y (input/output)
X = fft(x,N)*dt;
Y = fft(y,N)*dt;

% Calculate double sided CSD
% Will contain phase data
Gxy = (1/T)*(conj(X).*Y); % calculates Gxy (input/output)
Gyx = (1/T)*(conj(Y).*X); % calculates Gyx (output/input)

% Next steps turns Gxx into one sided CSD
Gxy = Gxy(1:ceil(N/2),:)*2;
Gyx = Gyx(1:ceil(N/2),:)*2;
for i = 1:Col(1,2)
    Gxy(1, i) = Gxy(1, i)/2;
    Gyx(1, i) = Gyx(1, i)/2;
end

% Average Gxy and Gyx (if multiple columns)
GxyM=mean(Gxy');
GyxM=mean(Gyx');

% Next determine frequency range
Fend=df*(N/2-1);
f=(0:df:Fend);

% plot Gxy and Gyx as subplots
subplot(211)
semilogy(f, abs(Gxy));
grid
title('Gxy Cross Spectrum Density, not averaged');
xlabel('Frequency, Hz');
ylabel('Cross Spectrum Density Mag.')

subplot(212)
semilogy(f, abs(Gyx));
grid
title('Gyx Cross Spectrum Density, not averaged');
xlabel('Frequency, Hz');
ylabel('Cross Spectrum Density Mag.')

% figure
% 
% subplot(211)
% semilogy(f, abs(GxyM));
% grid
% title('Gxy Cross Spectrum Density, averaged');
% xlabel('Frequency, Hz');
% ylabel('Cross Spectrum Density Mag.')
% 
% subplot(212)
% semilogy(f, abs(GyxM));
% grid
% title('Gyx Cross Spectrum Density, averaged');
% xlabel('Frequency, Hz');
% ylabel('Cross Spectrum Density Mag.')

figure

semilogy(f, abs(GxyM));
hold on
semilogy(f, abs(GyxM));
grid on

GxyOut = Gxy;
GyxOut = Gyx;
return