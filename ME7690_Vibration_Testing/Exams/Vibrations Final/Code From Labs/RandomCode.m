% Options1={'no', 1};
% 
% for i = 1:5;
% 
%     % get Tf for channel 2===============================================
%     [~, Random.ch2.TfH1(:,i)] = tfest(Random.ch2.Raw(:,i), ...
%         Random.ch3.Raw(:,i), dt, [], Options1);
% %     [~, Data.ch(2).TfcH1(:,i)] = tfest(Data.ch(2).Cor(:,i), ...
% %         Data.ch(3).Cor(:,i), dt, [], Options1);        
% %     %
% %     [~, Data.ch(2).TfrH2(:,i)] = tfest(Data.ch(2).Raw(:,i), ...
% %         Data.ch(3).Raw(:,i), dt, [], Options2);
% %     [~, Data.ch(2).TfcH2(:,i)] = tfest(Data.ch(2).Cor(:,i), ...
% %         Data.ch(3).Cor(:,i), dt, [], Options2);     
% %     
% %     % get Tf for channel 4===============================================
%     [~, Random.ch4.TfH1(:,i)] = tfest(Random.ch4.Raw(:,i), ...
%         Random.ch3.Raw(:,i), dt, [], Options1);
% %     [~, Data.ch(4).TfcH1(:,i)] = tfest(Data.ch(4).Cor(:,i), ...
% %         Data.ch(3).Cor(:,i), dt, [], Options1);        
% %     %
% %     [~, Data.ch(4).TfrH2(:,i)] = tfest(Data.ch(4).Raw(:,i), ...
% %         Data.ch(3).Raw(:,i), dt, [], Options2);
% %     [~, Data.ch(4).TfcH2(:,i)] = tfest(Data.ch(4).Cor(:,i), ...
% %         Data.ch(3).Cor(:,i), dt, [], Options2);  

%end

% for l = 1:length(NUMS);
% %     Data=[];
%     for p = 1:5;
% %         %Channel 2 ~output~
%         Ran_Tf2H1(l,p) = Random.ch2.TfH1(Q(l),p);
% %         Data.ch(2).TfcH1(:,p) = S(l).Data.ch(2).TfcH1(Q(l),p);
% %         Data.ch(2).TfrH2(:,p) = S(l).Data.ch(2).TfrH1(Q(l),p);
% %         Data.ch(2).TfcH2(:,p) = S(l).Data.ch(2).TfcH1(Q(l),p);
% %         
% %         %Channel 4 ~output~
%         Ran_Tf4H1(l,p) = Random.ch4.TfH1(Q(l),p);
% %         Data.ch(4).TfcH1(:,p) = S(l).Data.ch(4).TfcH1(Q(l),p);
% %         Data.ch(4).TfrH2(:,p) = S(l).Data.ch(4).TfrH1(Q(l),p);
% %         Data.ch(4).TfcH2(:,p) = S(l).Data.ch(4).TfcH1(Q(l),p);   
%      end
% %     Ext(l).Data = Data;
% %     Ext(l).num = NUMS(l);
% %     
% %     Tf2rH1(l, :) = Ext(l).Data.ch(2).TfrH1;
% %     Tf2cH1(l, :) = Ext(l).Data.ch(2).TfcH1;
% %     Tf2rH2(l, :) = Ext(l).Data.ch(2).TfrH2;
% %     Tf2cH2(l, :) = Ext(l).Data.ch(2).TfcH2;
% %     
% %     Tf4rH1(l, :) = Ext(l).Data.ch(4).TfrH1;
% %     Tf4cH1(l, :) = Ext(l).Data.ch(4).TfcH1;
% %     Tf4rH2(l, :) = Ext(l).Data.ch(4).TfrH2;
% %     Tf4cH2(l, :) = Ext(l).Data.ch(4).TfcH2;    
% %     
% end
% 
% Tf2rH1(:, 6) = mean(Tf2rH1')';
% Tf2cH1(:, 6) = mean(Tf2cH1')';
% Tf2rH2(:, 6) = mean(Tf2rH2')';
% Tf2cH2(:, 6) = mean(Tf2cH2')';
% %
% Tf4rH1(:, 6) = mean(Tf4rH1')';
% Tf4cH1(:, 6) = mean(Tf4cH1')';
% Tf4rH2(:, 6) = mean(Tf4rH2')';
% Tf4cH2(:, 6) = mean(Tf4cH2')';
% 
% %Get Coherence
% COH2r = abs(Tf2rH1(:, 6)./Tf2rH2(:, 6));
% COH4r = abs(Tf4rH1(:, 6)./Tf4rH2(:, 6));

% figure
% plot(Freq, 20*log10(mean(abs(Ran_Tf2H1'))'), 'linewidth', 2)
% hold on
% plot(Freq, 20*log10(mean(abs(Tf2rH1'))'), 'linewidth', 2)
% legend('FRF From Rand. Data', 'FRF From Sine Sweep', 'location', 'SouthEast');

figure
subplot(211);
plot(Freq, 20*log10(mean(abs(Ran_Tf2H1'))'), 'linewidth', 2);
title('FRF Estimate Channel 2 Random vs Sine Sweep');
xlabel('Frequency, Hz');
ylabel('Mag (dB)');
hold on
plot(Freq, 20*log10(abs(Tf2rH1(:,6))), 'linewidth', 2);
legend('FRF From Rand. Data', 'FRF From Sine Sweep', 'location', 'SouthEast');
grid on

subplot(212);
phase_Ran_Tf2H1 = unwrap(angle(mean(Ran_Tf2H1')'))*(180/pi);
plot(Freq, phase_Ran_Tf2H1, 'linewidth', 2);
title('Phase Estimate Channel 2 Random vs Sine Sweep');
xlabel('Frequency, Hz');
ylabel('Angle Deg.');
hold on
phase_Tf2rH1 = unwrap(angle(Tf2rH1(:,6)))*(180/pi);
plot(Freq, phase_Tf2rH1, 'linewidth', 2);
legend('FRF From Rand. Data', 'FRF From Sine Sweep', 'location', 'SouthWest');
grid on

figure
subplot(211);
plot(Freq, 20*log10(mean(abs(Ran_Tf4H1'))'), 'linewidth', 2);
title('FRF Estimate Channel 4 Random vs Sine Sweep');
xlabel('Frequency, Hz');
ylabel('Mag (dB)');
hold on
plot(Freq, 20*log10(abs(Tf4rH1(:,6))), 'linewidth', 2);
legend('FRF From Rand. Data', 'FRF From Sine Sweep', 'location', 'SouthEast');
grid on

subplot(212);
phase_Ran_Tf4H1 = unwrap(angle(mean(Ran_Tf4H1')'))*(180/pi);
plot(Freq, phase_Ran_Tf4H1, 'linewidth', 2);
title('Phase Estimate Channel 4 Random vs Sine Sweep');
xlabel('Frequency, Hz');
ylabel('Angle Deg.');
hold on
phase_Tf4rH1 = unwrap(angle(Tf4rH1(:,6)))*(180/pi);
plot(Freq, phase_Tf4rH1, 'linewidth', 2);
legend('FRF From Rand. Data', 'FRF From Sine Sweep', 'location', 'SouthWest');
grid on