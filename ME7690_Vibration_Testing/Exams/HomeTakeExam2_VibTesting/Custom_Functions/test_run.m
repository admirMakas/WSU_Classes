%clearvars Freq1 Freq2 Txf1 Txf2 i j

for i=1:2
    options{2}=i;
    [Freq1, Txf1(:,i)] = tfest(Chan_2_Filtered, Chan_1_Input_Shiffted, dt, options);
    %hold(subplot(211), 'on');
end
for j=1:2
    options{2}=j;
    %[Freq2, Txf2(:,j)] = tfest(Time_chan_2, Time_chan_1, 7.8125e-4, options);
    tfest(Chan_2_Filtered, Chan_1_Input_Shiffted, dt, options);
    hold(subplot(211), 'on');
end
figure
asd(Chan_2_Filtered,dt,4096)
figure
asd(Chan_1_Input_Shiffted,dt,4096)
figure
crsd(Chan_2_Filtered, Chan_1_Input_Shiffted,dt)
figure
coh(Chan_2_Filtered, Chan_1_Input_Shiffted,dt)