%% Code to plot the FRFs, Phase diagrams and COH for the sine sweep data
% The sine sweep goes from 0 to 55 HZ in 5 Hz increments.
%
% The experiment was conducted with 3 accelerometers on the wing. However,
% the acceleration sensor that pertains to channel #1 gave bad readings.
% Therefore, for the analysis only channel 2 and 4 were used. Channel 3
% measured the force input via the load cell.
%
% Raw data plots are for signals that were not modified in any before the
% FRF calculations.
%
% Cor data is corrected by removing the vertical shift that is present in
% the raw data. Function "ManFilt" was used to filter the data.
%%
% Following 2 plots show the FRF and phase plots for channel 2 RAW and COR
% data sets. The FRF plots have both H1 and H2 estimates ploted. However,
% they lie on top of each other. Removing the drift from the data has no
% significant impact on the FRF estimates.
%
%%
%Plot FRF for Channel 2 Raw H1 and H2
figure
subplot(211);
plot(Freq, 20*log10(abs(Tf2rH1(:, 6))), 'r-', 'linewidth', 2);
title('Mean FRF Estimate Channel 2 H1 and H2 Raw');
xlabel('Frequency, Hz');
ylabel('Mag (dB)');
hold on
plot(Freq, 20*log10(abs(Tf2rH2(:, 6))), 'k-', 'linewidth', 2);
grid on
subplot(212);
phase_Tf2rH1 = unwrap(angle(Tf2rH1(:, 6)))*(180/pi);
plot(Freq, phase_Tf2rH1, 'r-', 'linewidth', 2);
title('Phase Estimate Channel 2 Raw');
xlabel('Frequency, Hz');
ylabel('Angle Deg.');
grid on
%
%Plot FRF for Channel 2 Cor H1 and H2
figure
subplot(211);
plot(Freq, 20*log10(abs(Tf2cH1(:, 6))), 'r-', 'linewidth', 2);
title('Mean FRF Estimate Channel 2 H1 and H2 Cor');
xlabel('Frequency, Hz');
ylabel('Mag (dB)');
hold on
plot(Freq, 20*log10(abs(Tf2cH2(:, 6))), 'k-', 'linewidth', 2);
grid on
subplot(212);
phase_Tf2cH1 = unwrap(angle(Tf2cH1(:, 6)))*(180/pi);
plot(Freq, phase_Tf2cH1, 'r-', 'linewidth', 2);
title('Phase Estimate Channel 2 Cor');
xlabel('Frequency, Hz');
ylabel('Angle Deg.');
grid on
%
%==========================================================================
%%
% Next two plots show the FRF estimates for channel 4 data. It follows
% same trends that were observed for channel 2 data. Another thing that is
% obvious but should be mentioned is the fact that the sine sweep is not
% nearly sufficient to capture the response adequately. In reality the
% sweep should be far more refined. Also, there is about 3 points that
% define the peaks, which is not sufficient and can be seen in the phase
% plots as well.
%%
%Plot FRF for Channel 4 Raw H1 and H2
figure
subplot(211);
plot(Freq, 20*log10(abs(Tf4rH1(:, 6))), 'r-', 'linewidth', 2);
title('Mean FRF Estimate Channel 4 H1 and H2 Raw');
xlabel('Frequency, Hz');
ylabel('Mag (dB)');
hold on
plot(Freq, 20*log10(abs(Tf4rH2(:, 6))), 'k-', 'linewidth', 2);
grid on
subplot(212);
phase_Tf4rH1 = unwrap(angle(Tf4rH1(:, 6)))*(180/pi);
plot(Freq, phase_Tf4rH1, 'r-', 'linewidth', 2);
title('Phase Estimate Channel 4 Raw');
xlabel('Frequency, Hz');
ylabel('Angle Deg.');
grid on
%
%
%Plot FRF for Channel 4 Cor H1 and H2
figure
subplot(211);
plot(Freq, 20*log10(abs(Tf4cH1(:, 6))), 'r-', 'linewidth', 2);
title('Mean FRF Estimate Channel 4 H1 and H2 Cor');
xlabel('Frequency, Hz');
ylabel('Mag (dB)');
hold on
plot(Freq, 20*log10(abs(Tf4cH2(:, 6))), 'k-', 'linewidth', 2);
grid on
subplot(212);
phase_Tf4cH1 = unwrap(angle(Tf4cH1(:, 6)))*(180/pi);
plot(Freq, phase_Tf4cH1, 'r-', 'linewidth', 2);
title('Phase Estimate Channel 4 Cor');
xlabel('Frequency, Hz');
ylabel('Angle Deg.');
grid on
%
%==========================================================================
%%
% This section of plots shows channel 2 and 4 FRF estimates before all the
% data sets were averaged. 
%%
figure
plot(Freq, 20*log10(abs(Tf2rH1)), 'linewidth', 1);
title('FRF Estimates Channel 2 Raw');
xlabel('Frequency, Hz');
ylabel('Mag (dB)');
grid on
%
figure
plot(Freq, 20*log10(abs(Tf2cH1)), 'linewidth', 1);
title('FRF Estimates Channel 2 Cor');
xlabel('Frequency, Hz');
ylabel('Mag (dB)');
grid on
%
figure
plot(Freq, 20*log10(abs(Tf4rH1)), 'linewidth', 1);
title('FRF Estimates Channel 4 Raw');
xlabel('Frequency, Hz');
ylabel('Mag (dB)');
grid on
%
figure
plot(Freq, 20*log10(abs(Tf4cH1)), 'linewidth', 1);
title('FRF Estimates Channel 4 Cor');
xlabel('Frequency, Hz');
ylabel('Mag (dB)');
grid on
%
%==========================================================================
%%
% This section plots the coherence plots for channels 2 and 4.
%%
figure
plot(Freq, COH2r)
title('COH for Channel 2');
xlabel('Frequency, Hz');
ylabel('COH');
grid on
%
figure
plot(Freq, COH4r)
title('COH for Channel 2');
xlabel('Frequency, Hz');
ylabel('COH');
grid on