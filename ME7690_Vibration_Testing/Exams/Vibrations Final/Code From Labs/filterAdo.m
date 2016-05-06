%settings for beam data (channel 2) are the following:
%[0.01, 0.05, 60, 0.5]
d = fdesign.highpass('Fst,Fp,Ast,Ap', 0.01, 0.05, 30, 0.5);

designmethods(d)

Hd = design(d, 'equiripple');
fvtool(Hd);

n=Time_domain;
x=S(1).Data.ch(4).Raw(:,1);
y=filter(Hd,x);
Domega = (2*pi)/length(n);
freq = 0:(2*pi)/length(n):pi;
xdft = fft(x);
ydft = fft(y);
plot(freq, abs(xdft(1:length(x)/2+1)));
hold on;
plot(freq,abs(ydft(1:length(y)/2+1)), 'r', 'linewidth', 2);
legend('Original Signal', 'Filtered Signal', ...
    'Location', 'NorthEast');
title('Frequency Plot');
ylabel('Magnitude'); xlabel('Radians/Sample');
figure
plot(n, x, n, y);
title('Time History');
ylabel('Acceleration'); xlabel('Time (Sec)');