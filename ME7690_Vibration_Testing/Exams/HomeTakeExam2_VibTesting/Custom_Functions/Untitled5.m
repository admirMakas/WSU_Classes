load sensorData

t1 = (0:length(s1)-1)/Fs;
t2 = (0:length(s2)-1)/Fs;

subplot(2,1,1)
plot(t1,s1)
title('s_1')

subplot(2,1,2)
plot(t2,s2)
title('s_2')
xlabel('Time (s)')

[acor,lag] = xcorr(s2,s1);

[~,I] = max(abs(acor))
lagDiff = lag(I)
timeDiff = lagDiff/Fs

figure
plot(lag,acor)
% a3 = gca;
% a3.XTick = sort([-3000:1000:3000 lagDiff]);