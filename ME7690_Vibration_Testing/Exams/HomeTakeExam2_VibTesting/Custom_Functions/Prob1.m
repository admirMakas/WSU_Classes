%%
% Use crcor function to get sigal correlations
[t1, yout_xy] = crcor(xn, yn, dt, 1);
[t1, yout_yx] = crcor(yn, xn, dt, 1);

%%
% %Plot time history
% plot(t, xn)
% axis([0 1.05 -inf inf]);
% legend('X Data', 'Location', 'Southeast');
% title('Time History');
% xlabel('Time, Sec');
% ylabel('Mag');
% 
% figure
% plot(t, yn)
% axis([0 1.05 -inf inf]);
% legend('Y Data', 'Location', 'Southeast');
% title('Time History');
% xlabel('Time, Sec');
% ylabel('Mag');

%Plot cross correlation for XY and YX
figure
plot(t1, yout_xy, t1, yout_yx)
axis([-1.1 1.1 -inf inf]);
legend('XY Corr', 'YX Corr', 'Location', 'Southeast');
title('Signal Correlations');
xlabel('Time, Sec');
ylabel('Mag');

%%
%Plot both signals to obtain natural frequencies
FFTPlot(xn, dt);
FFTPlot(yn, dt);

%%
% Use matlab function xcorr to plot cross correlation
[acor, lag] = xcorr(yn, xn);

[~, I] = max(abs(acor))
lagDiff = lag(I)
timeDiff = lagDiff/1025

%plot(lag,acor)