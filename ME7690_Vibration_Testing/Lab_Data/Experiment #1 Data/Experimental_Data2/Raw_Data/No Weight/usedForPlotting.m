plot(Time,loc4_output_chan_wt_weight_raw(:,1));
hold on
plot(Time,loc4_output_chan_wt_weight_corr(:,1));
legend('Raw Data', 'Corrected Data')
title('Time History for LOC4');
xlabel('Time (Sec)');
ylabel('Acceleration');
grid on
%
%
figure
plot(Time,base_output_chan_wt_weight_raw(:,1));
hold on
plot(Time,base_output_chan_wt_weight_corr(:,1));
legend('Raw Data', 'Corrected Data')
title('Time History for Base');
xlabel('Time (Sec)');
ylabel('Acceleration');
grid on
%
%
figure
plot(Time,base_input_chan_wt_weight_raw(:,1));
hold on
plot(Time,base_input_chan_wt_weight_corr(:,1));
legend('Raw Data', 'Corrected Data')
title('Time History for Base Input');
xlabel('Time (Sec)');
ylabel('Magnitude');
grid on
%
%
figure
plot(Time,loc4_input_chan_wt_weight_raw(:,1));
hold on
plot(Time,loc4_input_chan_wt_weight_corr(:,1));
legend('Raw Data', 'Corrected Data')
title('Time History for LOC4 Input');
xlabel('Time (Sec)');
ylabel('Magnitude');
grid on