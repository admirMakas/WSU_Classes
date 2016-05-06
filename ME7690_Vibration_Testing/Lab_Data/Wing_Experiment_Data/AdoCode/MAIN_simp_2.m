%% Lab #5                                      % Student:      Daniel Clark
clear all                                     % Instructor: Dr. Ha-Rok Bae
close all                                      % Class:             ME-4210
clc
format short
fontSize = 16.0;

% Noise on signals 
 


for i = 1:5
    FileName = strcat('NoInput_',int2str(i),'.mat');
    load(FileName,'-mat')
    
    Noise.ch1.Raw(:,i) = Time_chan_1;
    Noise.ch2.Raw(:,i) = Time_chan_2;
    Noise.ch3.Raw(:,i) = Time_chan_3;
    Noise.ch4.Raw(:,i) = Time_chan_4;
    
end

Noise.ch1.Mean = mean(Noise.ch1.Raw')'; Noise.ch2.Mean  = mean(Noise.ch2.Raw')';
Noise.ch3.Mean  = mean(Noise.ch3.Raw')'; Noise.ch4.Mean  = mean(Noise.ch4.Raw')';
Noise.t = Time_domain;
Noise.f = Freq_domain;

figure, hold on
plot(Noise.t, Noise.ch1.Raw(:,1), Noise.t, Noise.ch1.Raw(:,2), Noise.t, Noise.ch1.Raw(:,3), ...
    Noise.t, Noise.ch1.Raw(:,4), Noise.t, Noise.ch1.Mean)
plot(Noise.t, Noise.ch1.Mean,'k','LineWidth',2)

figure, hold on
plot(Noise.t, Noise.ch1.Mean,'k','LineWidth',2)
plot(Noise.t, Noise.ch2.Mean,'b','LineWidth',2)
plot(Noise.t, Noise.ch3.Mean,'r','LineWidth',2)
plot(Noise.t, Noise.ch4.Mean,'--k','LineWidth',2)


for i = 1:5
    %FileName = strcat('Random_',int2str(i),'.mat');
    FileName = strcat('freq_0025_00_',int2str(i),'.mat');
    load(FileName,'-mat')
    
    Data.ch(1).Raw(:,i) = Time_chan_1;
    Data.ch(1).Cor(:,i) = Time_chan_1 - Noise.ch1.Mean;
    Data.ch(2).Raw(:,i) = Time_chan_2;
    Data.ch(2).Cor(:,i) = Time_chan_2 - Noise.ch2.Mean;
    Data.ch(3).Raw(:,i) = Time_chan_3;
    Data.ch(3).Cor(:,i) = Time_chan_3 - Noise.ch3.Mean;
    Data.ch(4).Raw(:,i) = Time_chan_4;
    Data.ch(4).Cor(:,i) = Time_chan_4 - Noise.ch4.Mean;
    
end

[ S ] = SignalAnalysis( Data, Noise.t);
figure, hold on
plot(S.f, S.Sensor(4).H(1).mag)
plot(S.f, S.Sensor(4).H(2).mag)
plot(S.f, S.Sensor(4).H(3).mag)



figure, hold on
plot(S.f, S.Sensor(1).FFT)
plot(S.f, S.Sensor(2).FFT)
plot(S.f, S.Sensor(4).FFT)

figure, hold on
plot(S.f, S.Sensor(3).FFT)

break

% plot(S.f, S.Sensor(2).H(2).mag)
% plot(S.f, S.Sensor(3).H(2).mag,'k','LineWidth',2)
% plot(S.f, S.Sensor(4).H(2).mag)

figure, hold on
plot(S.f, S.Sensor(1).H(1).phase)
plot(S.f, S.Sensor(2).H(2).phase)
plot(S.f, S.Sensor(3).H(3).phase,'k','LineWidth',2)
%plot(S.f, S.Sensor(4).H(1).phase)

figure, hold on
plot(S.f, S.Sensor(1).Coh)
plot(S.f, S.Sensor(2).Coh)
plot(S.f, S.Sensor(3).Coh,'k','LineWidth',2)
plot(S.f, S.Sensor(4).Coh)


break
MEAN = mean(Data.ch1.Raw')';


figure, hold on
plot(Noise.t, Data.ch(1).Raw(:,1), Noise.t, Data.ch(1).Raw(:,2), Noise.t, Data.ch(1).Raw(:,3), ...
    Noise.t, Data.ch(1).Raw(:,4),Noise.t, MEAN,'k','linewidth',2)
%plot(Noise.t, Data.ch1.Mean,'k','LineWidth',2)

break


load('NoInput_1.mat','-mat')
noise = Time_chan_1;
noiseBeam = Time_chan_2;
noiseBase = Time_chan_4;
figure
plot(Time_domain, Time_chan_4)



% initialization 
i = 1;
FileName = strcat('Case_base_W_',int2str(i),'.mat');
load(FileName,'-mat')

% Base
out_base_noW_accMatrix_currected = zeros(length(Time_chan_1),5);
out_base_noW_accMatrix = out_base_noW_accMatrix_currected;
in_base_noW_accMatrix_currected = out_base_noW_accMatrix_currected;
in_base_noW_accMatrix = out_base_noW_accMatrix_currected;

% position 4
out_pos4_noW_accMatrix_currected = out_base_noW_accMatrix_currected;
out_pos4_noW_accMatrix = out_base_noW_accMatrix_currected;
in_pos4_noW_accMatrix_currected = out_base_noW_accMatrix_currected;
in_pos4_noW_accMatrix = out_base_noW_accMatrix_currected;

% load the rest of the data
for i = 1:5
    FileName = strcat('Case_base_W_',int2str(i),'.mat');
    load(FileName,'-mat')
    out_base_noW_accMatrix_currected(:,i) = Time_chan_2 - noiseBase;
    out_base_noW_accMatrix(:,i) = Time_chan_2;
    in_base_noW_accMatrix_currected(:,i) = Time_chan_1 - noiseHammer;
    in_base_noW_accMatrix(:,i) = Time_chan_1;
    
    FileName = strcat('Case_pos4_noW_',int2str(i),'.mat');
    load(FileName,'-mat')
    out_pos4_noW_accMatrix_currected(:,i) = Time_chan_2 - noiseBase;
    out_pos4_noW_accMatrix(:,i) = Time_chan_2;
    in_pos4_noW_accMatrix_currected(:,i) = Time_chan_1 - noiseHammer;
    in_pos4_noW_accMatrix(:,i) = Time_chan_1;
end  
% Averaging
out_base_noW_accMatrix_currected_avg = mean(out_base_noW_accMatrix_currected')';
out_base_noW_accMatrix_avg = mean(out_base_noW_accMatrix')';
in_base_noW_accMatrix_currected_avg = mean(in_base_noW_accMatrix_currected')';
in_base_noW_accMatrix_avg = mean(in_base_noW_accMatrix')';

out_pos4_noW_accMatrix_currected_avg = mean(out_pos4_noW_accMatrix_currected')';
out_pos4_noW_accMatrix_avg = mean(out_pos4_noW_accMatrix')';
in_pos4_noW_accMatrix_currected_avg = mean(in_pos4_noW_accMatrix_currected')';
in_pos4_noW_accMatrix_avg = mean(in_pos4_noW_accMatrix')';


subplot(2,1,1)
hold on
plot(Time_domain,out_base_noW_accMatrix_avg)
plot(Time_domain,noiseBeam)
grid
axis([0 3 min(out_base_noW_accMatrix_avg) max(out_base_noW_accMatrix_avg)])
xlabel(' Time (sec) ')
ylabel(' Acceleration (m/s^2)')

subplot(2,1,2) 
hold on
plot(Time_domain,out_base_noW_accMatrix_currected_avg)
grid
axis([0 3 min(out_base_noW_accMatrix_currected_avg) max(out_base_noW_accMatrix_currected_avg)])
xlabel(' Time (sec) ')
ylabel(' Acceleration (m/s^2) ')

figure
hold on
plot(Time_domain,out_base_noW_accMatrix_avg)
plot(Time_domain,out_base_noW_accMatrix_currected_avg)
grid
axis([0 3 min(out_base_noW_accMatrix_avg) max(out_base_noW_accMatrix_avg)])
axis([0 3 -70 70])
xlabel(' Time (sec) ','fontSize',fontSize);
ylabel(' Acceleration (m/s^2) ','fontSize',fontSize);
set(gca,'fontSize',fontSize);

figure
hold on
plot(Time_domain,out_pos4_noW_accMatrix_avg)
plot(Time_domain,out_pos4_noW_accMatrix_currected_avg)
grid
axis([0 3 -15 15])
xlabel(' Time (sec) ','fontSize',fontSize);
ylabel(' Acceleration (m/s^2) ','fontSize',fontSize);
set(gca,'fontSize',fontSize);


t = Time_domain;
w = Freq_domain;
x = out_base_noW_accMatrix;
f = in_base_noW_accMatrix;

[freqout,H1]=tfest(x,f,t);
OPTIONS{1} = 'yes';
OPTIONS{2} = 2;
[freqout,H2]=tfest(x,f,t,[],OPTIONS);
OPTIONS{1} = 'yes';
OPTIONS{2} = 3;
[freqout,Hv]=tfest(x,f,t,[],OPTIONS);

figure
hold on
plot(freqout,20*log10(abs(H1)),'r','LineWidth',2)
plot(freqout,20*log10(abs(H2)),'b','LineWidth',2)
plot(freqout,20*log10(abs(Hv)),'k','LineWidth',2)
xlabel(' Frequency (HZ) ','fontSize',fontSize);
ylabel(' Magnitude (dB) ','fontSize',fontSize);
set(gca,'fontSize',fontSize);
grid on

t = Time_domain;
w = Freq_domain;
x = out_base_noW_accMatrix_currected;
f = in_base_noW_accMatrix_currected;

[H1_mine, H2_mine, wnew,H1Matrix,H2Matrix] = Find_H1_H2( t, w, x, f );
%plot(wnew,20*log10(abs(H1)),'g')
figure
hold on
plot(wnew,H1_mine.mag,'g')
plot(wnew,H2_mine.mag,'k')
plot(freqout,20*log10(abs(H1)),'r','LineWidth',2)
plot(freqout,20*log10(abs(H2)),'b','LineWidth',2)

figure
hold on
plot(wnew,H1_mine.phase,'g')
plot(wnew,H2_mine.phase,'k')
plot(wnew,unwrap(angle(H1)) * 180 / pi,'r')
plot(wnew,unwrap(angle(H2)) * 180 / pi,'b')


figure
hold on
plot(wnew,H1_mine.phase,'r')
plot(wnew,H2_mine.phase,'b')
% hold on
% [Gxx] = AutoSD(x(:,5),t(2))
break
clearvars -except H1_mine H2_mine H1 H2

break
%% Performing my FFT

t = Time_domain;
w = Freq_domain;
x = out_base_noW_accMatrix_currected_avg;
f = in_base_noW_accMatrix_currected_avg;

% clearvars -except t w x f

dt = t(end) / length(t);
dw = w(end) / length(w);


%for i = 1:5


F = fft(f) * dt;
F = F(1:length(w));
X = fft(x) * dt;
X = X(1:length(w));


Gfx = 2 * (conj(F) .* X) * dw; Gfx(1) = Gfx(1) / 2;
Gxf = 2 * (conj(X) .* F) * dw; Gxf(1) = Gxf(1) / 2;
Gff = 2 * real(conj(F) .* F) * dw; Gff(1) = Gff(1) / 2;
Gxx = 2 * real(conj(X) .* X) * dw; Gxx(1) = Gxx(1) / 2;

H1 = Gfx ./ Gff;
H2 = Gxx ./ Gxf;

mag_H1 = 20 * log10(abs(H1));
phase_H1 = angle(H1) * 180 / pi;
mag_H2 = 20 * log10(abs(H2));
phase_H2 = angle(H2) * 180 / pi;

t = Time_domain;
w = Freq_domain;
x = out_base_noW_accMatrix_currected;
f = in_base_noW_accMatrix_currected;

[H1, H2, wnew] = Find_H1_H2( t, w, x, f );

figure
hold on
plot(wnew,H1.mag,'r')
plot(wnew,H2.mag,'b')
plot(wnew,mag_H1,'g')
plot(wnew,mag_H2,'k')

figure
subplot(2,1,1)
hold on
plot(w,H1.mag)
plot(w,H2.mag)
grid
%axis([0 3 min(out_base_noW_accMatrix_avg) max(out_base_noW_accMatrix_avg)])
xlabel(' Time (sec) ')
ylabel(' Acceleration (m/s^2)')

subplot(2,1,2) 
hold on
plot(w,H1.phase)
plot(w,H2.phase)
grid
%axis([0 3 min(out_base_noW_accMatrix_currected_avg) max(out_base_noW_accMatrix_currected_avg)])
xlabel(' Time (sec) ')
ylabel(' Acceleration (m/s^2) ')




break
H = 20 * log10(abs(Hf_chan_2));
H1bob = conj(Hf_Cross_Spec_chan_2) ./ PSD_chan_1; % From BOBCAT
H2bob = PSD_chan_2 ./ Hf_Cross_Spec_chan_2; % From BOBCAT


break





%% Using Metric units
% Givens
L = 21.75*0.0254;            
H = 0.5*0.0254;             
W = 1*0.0254;              
E = 7.1*10^10;            

% Beam Analytical
display('*** Beam Analytical ***')
density = 2700;        
volume = L*H*W;         
I = (1/12)*W*H^3;
m = density * volume    
k = (3*E*I) / (L^3)
w_n = sqrt(k/m)        

% experimental
display('*** Experimental ***')
load case1.mat
hold on
acceleration = Time_chan_2 - Time_chan_4;

time = Time_domain(length(Time_domain));
n = 0;
for i = 2:length(acceleration)
    if sign(acceleration(i)) == sign(acceleration(i-1))
        
    else
       n = n +1;
    end
end


n = 107; % number of peaks
T_d = time/n
w_d = 2*pi/T_d
position = acceleration / (-(w_d^2));
log_dec = -(1/n)*log((9.311*10^-6)/(8.719*10^-5)) % where do these numbers come from...
damp_ratio = log_dec / (sqrt(4*pi^2 + log_dec^2))
W_n = w_d / sqrt(1-damp_ratio^2)
cr = 2 * sqrt(k*m);
c_exp = cr*damp_ratio

% Analytical again 
display('*** Beam Analytical ***')
c_analytical = w_n*damp_ratio

% plot data
subplot(2,1,1)
hold on
plot(Time_domain,acceleration)
acceleration3 = Time_chan_4;
plot(Time_domain,acceleration3)
grid
axis([0 3 min(acceleration) max(acceleration)])
title(' Case 1 Data ' )
xlabel(' Time (sec) ')
ylabel(' Acceleration (m/s^2)')

subplot(2,1,2) 
hold on
plot(Time_domain,position)
grid
axis([0 3 min(position) max(position)])
xlabel(' Time (sec) ')
ylabel(' Displacement (m)')
plot(Time_domain,max(position).*exp(Time_domain*-log_dec*250),'k','LineWidth',2)

fontSize = 18.0;
lineWidth = 6.0;
% -----------------
x = Time_chan_2;
f = Time_chan_1;
t = Time_domain;
w = Freq_domain;
dt = t(end) / length(t);
dw = w(end) / length(w);
% ----------------------------------------------------------------------- %
F = fft(f) * dt; F = F(1:length(w));
X = fft(x) * dt; X = X(1:length(w));


Gfx = 2 * (conj(F) .* X) * dw; Gfx(1) = Gfx(1) / 2;
Gxf = 2 * (conj(X) .* F) * dw; Gxf(1) = Gxf(1) / 2;
Gff = 2 * real(conj(F) .* F) * dw; Gff(1) = Gff(1) / 2;
Gxx = 2 * real(conj(X) .* X) * dw; Gxx(1) = Gxx(1) / 2;

H1 = Gfx ./ Gff;
H2 = Gxx ./ Gxf;
H = 20 * log10(abs(Hf_chan_2));
H1bob = conj(Hf_Cross_Spec_chan_2) ./ PSD_chan_1; % From BOBCAT
H2bob = PSD_chan_2 ./ Hf_Cross_Spec_chan_2; % From BOBCAT

[H1Mag H1Phase] = FRFpost(H1);
[H2Mag H2Phase] = FRFpost(H2);
[H1bobMag H1bobPhase] = FRFpost(H1bob);
[H2bobMag H2bobPhase] = FRFpost(H2bob);
[HMag HPhase] = FRFpost(Hf_chan_2);

figure
FRFplot(w(1),w(end),H1bob)

figure
FRFplot(w(1),w(end),H2bob)


figure
% plot(w,H1Mag(1:length(w)),'r',w,H2Mag(1:length(w)),'k','lineWidth',lineWidth)
plot(w,H1Mag,'r',w,H2Mag,'k','lineWidth',lineWidth)
title('H_1 and H_2 from G_{xx}, G_{ff}, G_{xf}, G_{fx}')
legend('H_1','H_2')
xlabel('Freq, Hz','fontSize',fontSize);
ylabel('Mag, dB','fontSize',fontSize);
set(gca,'fontSize',fontSize);

figure
plot(w,H1bobMag,'r',w,H2bobMag,'k','lineWidth',lineWidth)
legend('H_1','H_2')
xlabel('Freq, Hz','fontSize',fontSize);
ylabel('Mag, dB','fontSize',fontSize);
set(gca,'fontSize',fontSize);

figure
plot(w,Hf_coh_chan_2,'k','lineWidth',lineWidth)
xlabel('Freq, Hz','fontSize',fontSize);
ylabel('Coherence','fontSize',fontSize);
set(gca,'fontSize',fontSize);
