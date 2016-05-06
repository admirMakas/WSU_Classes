function [Gxx] = AutoSD(x,dt)

% get data size
N = length(x);

% determine # of columns
Col = size(x);

% calculate time period
T = N*dt;

% calculate frequency step
df = 1/T;

% put data into frequency domain
X = fft(x,N)*dt;

% Calculate double sided ASD
% Does not contain phase data
Gxx = (1/T)*real(conj(X).*X);

% Next steps turns Gxx into one sided ASD
Gxx = Gxx(1:ceil(N/2),:)*2;
for i = 1:Col(1,2)
    Gxx(1, i) = Gxx(1, i)/2;
end

% Average Gxx (if multiple columns)
GxxM=mean(Gxx')';

% Next determine frequency range
Fend=df*(N/2-1);
f=(0:df:Fend)';

% plot test
if Col(1,2)~=1
    plot(f, 20*log10(GxxM));
end
figure
plot(f, 20*log10(Gxx))

GxxOut = Gxx;
return