%% Lab #5                                      % Student:      Daniel Clark
clear all                                     % Instructor: Dr. Ha-Rok Bae
close all                                      % Class:             ME-4210
clc
format short
fontSize = 16.0;

% Noise on signals 

for i = 1:5;
    FileName = strcat('NoInput_',int2str(i),'.mat');
    load(FileName,'-mat');
    
    Noise.ch1.Raw(:,i) = Time_chan_1;
    Noise.ch2.Raw(:,i) = Time_chan_2;
    Noise.ch3.Raw(:,i) = Time_chan_3;
    Noise.ch4.Raw(:,i) = Time_chan_4;
    
end

Noise.ch1.Mean = mean(Noise.ch1.Raw')'; Noise.ch2.Mean  = mean(Noise.ch2.Raw')';
Noise.ch3.Mean = mean(Noise.ch3.Raw')'; Noise.ch4.Mean  = mean(Noise.ch4.Raw')';
Noise.t = Time_domain;
Noise.f = Freq_domain;
NUMS = 5:5:55;

S = [];
for j = 1:length(NUMS);
Data = [];
for i = 1:5;
    if NUMS(j) == 20 && i ==3;
        
    else
    FileName = strcat('freq_00',figureTag(NUMS(j)),'_00_',int2str(i),'.mat');
    load(FileName,'-mat');
    
    Data.ch(1).Raw(:,i) = Time_chan_1;
    Data.ch(1).Cor(:,i) = Time_chan_1 - Noise.ch1.Mean;
    Data.ch(2).Raw(:,i) = Time_chan_2;
    Data.ch(2).Cor(:,i) = Time_chan_2 - Noise.ch2.Mean;
    Data.ch(3).Raw(:,i) = Time_chan_3;
    Data.ch(3).Cor(:,i) = Time_chan_3 - Noise.ch3.Mean;
    Data.ch(4).Raw(:,i) = Time_chan_4;
    Data.ch(4).Cor(:,i) = Time_chan_4 - Noise.ch4.Mean;
    end
end

S(j).Data = Data;
%S(j).obj  = SignalAnalysis( Data, Noise.t);
S(j).num = NUMS(j);

end




