% x = input
% y = output
function [TfuncH1, TfuncH2] = Tfunc(x, y, dt)

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

% Calculate double sided ASD
% no phase data
Gxx = (1/T)*real(conj(X).*X); % input ASD
Gyy = (1/T)*real(conj(Y).*Y); % output ASD

% Calculate double sided CSD
% Will contain phase data
Gxy = (1/T)*(conj(X).*Y); % calculates Gxy (input/output)
Gyx = (1/T)*(conj(Y).*X); % calculates Gyx (output/input)

% Next steps turns Gxx into one sided ASD
Gxx = Gxx(1:ceil(N/2),:)*2;
Gyy = Gyy(1:ceil(N/2),:)*2;
for i = 1:Col(1,2)
    Gxx(1, i) = Gxx(1, i)/2;
    Gyy(1, i) = Gyy(1, i)/2;
end

GxxM=mean(Gxx');
GyyM=mean(Gyy');

% Next steps turns Gxx into one sided CSD
Gxy = Gxy(1:ceil(N/2),:)*2;
Gyx = Gyx(1:ceil(N/2),:)*2;
for i = 1:Col(1,2)
    Gxy(1, i) = Gxy(1, i)/2;
    Gyx(1, i) = Gyx(1, i)/2;
end

GxyM=mean(Gxy');
GyxM=mean(Gyx');

% Next calculate H1 estimate
TfuncH1 = Gyx./Gxx;
TfuncH1M = mean(TfuncH1');
%TfuncH1M = conj(TfuncH1M);

% Next calculate H2 estimate
TfuncH2 = Gyy./Gyx;
TfuncH2M = mean(TfuncH2');
%TfuncH2M = conj(TfuncH2M);

% Next determine frequency range
Fend=df*(N/2-1);
f=(0:df:Fend);

% plot H1 and phase as subplots
subplot(211)
plot(f, 20*log10(abs(TfuncH1M)));
grid
title('H1 estimate, averaged');
xlabel('Frequency, Hz');
ylabel('H1 Mag.')
% plot phase
subplot(212)
phaseH1 = unwrap(angle(TfuncH1M))*(180/pi);
plot(f, phaseH1);
grid

figure

% plot H2 and phase as subplots
subplot(211)
plot(f, 20*log10(abs(TfuncH2M)));
grid
title('H2 estimate, averaged');
xlabel('Frequency, Hz');
ylabel('H2 Mag.')
% plot phase
subplot(212)
phaseH2 = unwrap(angle(TfuncH2M))*(180/pi);
plot(f, phaseH2);
grid

figure

plot(f, 20*log10(abs(TfuncH1M)));
hold on
plot(f, 20*log10(abs(TfuncH2M)));
grid on

H1Out = TfuncH1M;
H2Out = TfuncH2M;
return