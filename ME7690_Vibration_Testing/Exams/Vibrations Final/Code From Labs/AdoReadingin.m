%% Lab #5                                      % Student:      Daniel Clark
clear all                                     % Instructor: Dr. Ha-Rok Bae
close all                                      % Class:             ME-4210
clc
format short
fontSize = 16.0;

%==========================================================================
%Define manual filter that will in this case remove the first 5 Fourier
%coefficients in the frequency domain and use ifft to get a new time series
%that will remove some of the noise. This new dataset will be called "Cor"
%==========================================================================

%Right now must be defined in the code!!!!
%used for function ManFilt call.
df = 0.3125;
dt=0.00039062501;
ave='noave';
Freq = [5 10 13.50 15 20 25 30 35 40 45 50 51.70 55];

% Random data signals 

for i = 1:5;
    FileName = strcat('Random_',int2str(i),'.mat');
    load(FileName,'-mat');
    
    Random.ch1.Raw(:,i) = Time_chan_1;
    Random.ch2.Raw(:,i) = Time_chan_2;
    Random.ch3.Raw(:,i) = Time_chan_3;
    Random.ch4.Raw(:,i) = Time_chan_4;
    
end

% Noise.ch1.Mean = mean(Noise.ch1.Raw')'; Noise.ch2.Mean  = mean(Noise.ch2.Raw')';
% Noise.ch3.Mean = mean(Noise.ch3.Raw')'; Noise.ch4.Mean  = mean(Noise.ch4.Raw')';
% Noise.t = Time_domain;
% Noise.f = Freq_domain;

TimeRange = Time_domain;
FreqRange = Freq_domain;

NUMS = 5:5:55;
NUMS = [NUMS 13 51];
NUMS=sort(NUMS);

NUMS2 = {'00', '00', '50', '00', '00', '00', '00', '00', '00', ...
    '00', '00', '70', '00'}';

S = [];
for j = 1:length(NUMS);
Data = [];
for i = 1:5;
    if NUMS(j) == 20 && i ==3;
        
    else
    FileName = strcat('freq_00',figureTag(NUMS(j)),'_',NUMS2(j,1),'_',int2str(i),'.mat');
    FileName=char(FileName{1})
    load(FileName,'-mat');
    
    Data.ch(1).Raw(:,i) = Time_chan_1;
    Data.ch(1).Cor(:,i) = ManFilt(Time_chan_1, dt);
    %
    Data.ch(2).Raw(:,i) = Time_chan_2;
    Data.ch(2).Cor(:,i) = ManFilt(Time_chan_2, dt);
    %
    % get ASD==============================================================
    [~, Data.ch(2).Gffr(:,i)] = asd(Data.ch(2).Raw(:,i), dt, [], ave);
    [~, Data.ch(2).Gffc(:,i)] = asd(Data.ch(2).Cor(:,i), dt, [], ave);
    FreqRange = asd(Data.ch(2).Cor(:,i), dt, [], ave);
    % =====================================================================
    %
    Data.ch(3).Raw(:,i) = Time_chan_3;
    Data.ch(3).Cor(:,i) = ManFilt(Time_chan_3, dt);
    %
    % get ASD==============================================================
    [~, Data.ch(3).Gffr(:,i)] = asd(Data.ch(3).Raw(:,i), dt, [], ave);
    [~, Data.ch(3).Gffc(:,i)] = asd(Data.ch(3).Cor(:,i), dt, [], ave);    
    %======================================================================
    %
    Data.ch(4).Raw(:,i) = Time_chan_4;
    Data.ch(4).Cor(:,i) = ManFilt(Time_chan_4, dt);
    %
    % get ASD==============================================================
    [~, Data.ch(4).Gffr(:,i)] = asd(Data.ch(4).Raw(:,i), dt, [], ave);
    [~, Data.ch(4).Gffc(:,i)] = asd(Data.ch(4).Cor(:,i), dt, [], ave);    
    %======================================================================
    
    % get CRSD for channel 2===============================================
    [~, Data.ch(2).Gxfr(:,i)] = crsd(Data.ch(2).Raw(:,i), Data.ch(3).Raw(:,i), dt, [], ave);
    [~, Data.ch(2).Gxfc(:,i)] = crsd(Data.ch(2).Cor(:,i), Data.ch(3).Cor(:,i), dt, [], ave);        
    %
    [~, Data.ch(2).Gfxr(:,i)] = crsd(Data.ch(3).Raw(:,i), Data.ch(2).Raw(:,i), dt, [], ave);
    [~, Data.ch(2).Gfxc(:,i)] = crsd(Data.ch(3).Cor(:,i), Data.ch(2).Cor(:,i), dt, [], ave);        
    %======================================================================    
    
    % get CRSD for channel 4===============================================
    [~, Data.ch(4).Gxfr(:,i)] = crsd(Data.ch(4).Raw(:,i), Data.ch(3).Raw(:,i), dt, [], ave);
    [~, Data.ch(4).Gxfc(:,i)] = crsd(Data.ch(4).Cor(:,i), Data.ch(3).Cor(:,i), dt, [], ave);        
    %
    [~, Data.ch(4).Gfxr(:,i)] = crsd(Data.ch(3).Raw(:,i), Data.ch(4).Raw(:,i), dt, [], ave);
    [~, Data.ch(4).Gfxc(:,i)] = crsd(Data.ch(3).Cor(:,i), Data.ch(4).Cor(:,i), dt, [], ave);        
    %======================================================================        

    end
end

S(j).Data = Data;
%S(j).obj  = SignalAnalysis( Data, Noise.t);
S(j).num = NUMS(j);

end

%Find index of frequency values needed to constuct transfer function
Q=[];
for k = 1:length(Freq);
    B=(FreqRange>(Freq(k)-0.15) & FreqRange<(Freq(k)+0.15));
    Qi=find(B==1);
    Q = [Q; Qi];
end

Ext=[];
%Extract values from the PSD and CRSD vectors
for l = 1:length(NUMS);
    %Data=[];
    for p = 1:5;
        %Channel 2 ~output~
        Gxxr2(:,p) = S(l).Data.ch(2).Gffr(Q,p);
        Gxxc2(:,p) = S(l).Data.ch(2).Gffc(Q,p);
        %
        Gxfr2(:,p) = S(l).Data.ch(2).Gxfr(Q,p);
        Gxfc2(:,p) = S(l).Data.ch(2).Gxfc(Q,p);
        Gfxr2(:,p) = S(l).Data.ch(2).Gfxr(Q,p);
        Gfxc2(:,p) = S(l).Data.ch(2).Gfxc(Q,p);
        %
        %Channel 4 ~output~
        Gxxr4(:,p) = S(l).Data.ch(4).Gffr(Q,p);
        Gxxc4(:,p) = S(l).Data.ch(4).Gffc(Q,p);
        %
        Gxfr4(:,p) = S(l).Data.ch(4).Gxfr(Q,p);
        Gxfc4(:,p) = S(l).Data.ch(4).Gxfc(Q,p);
        Gfxr4(:,p) = S(l).Data.ch(4).Gfxr(Q,p);
        Gfxc4(:,p) = S(l).Data.ch(4).Gfxc(Q,p);   
        %
        %Channel 3 ~input~
        Gffr3(:,p) = S(l).Data.ch(3).Gffr(Q,p);
        Gffc3(:,p) = S(l).Data.ch(3).Gffc(Q,p);        
    end
    %Ext.Data = Data;
end

clearvars -except S ave dt ave Freq FreqRange Q Gxxr2 Gxxc2 Gxfr2 Gxfc2 ...
    Gfxr2 Gfxc2 Gxxr4 Gxxc4 Gxfr4 Gxfc4 Gfxr4 Gfxc4 Gffr3 Gffc3 TimeRange ...
    NUMS NUMS2 df Random

H1_2r=Gfxr2./Gffr3;
H1_2c=Gfxc2./Gffc3;

H2_2r = Gxxr2./Gxfr2;
H2_2c = Gxxc2./Gxfc2;

H1_4r=Gfxr4./Gffr3;
H1_4c=Gfxc4./Gffc3;

H2_4r = Gxxr4./Gxfr4;
H2_4c = Gxxc4./Gxfc4;

%plot(Freq', 20*log(abs(mean(H1_2r')')))

plot(Freq', 20*log(abs((H2_4c))))