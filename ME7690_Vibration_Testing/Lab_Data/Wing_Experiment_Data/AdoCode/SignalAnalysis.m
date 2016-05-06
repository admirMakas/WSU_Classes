function [ STRUCTURE ] = SignalAnalysis( Data, t)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


f = Data.ch(3).Cor;

freqout = [];

for i = 1:4
x = Data.ch(i).Raw;
X = fft(x(:,1),length(t)) * (t(2)-t(1));
X = X(1:length(X)/2);
X(1) = 0;


[~,H]=tfest(x,f,t);
STRUCTURE.Sensor(i).H(1).total = H;
[mag, phase ] = mag_phase(H);
STRUCTURE.Sensor(i).H(1).mag = mag;
STRUCTURE.Sensor(i).H(1).phase = phase;

OPTIONS{1} = 'yes';
OPTIONS{2} = 2;
[~,H]=tfest(x,f,t,[],OPTIONS);
STRUCTURE.Sensor(i).H(2).total = H;
[mag, phase ] = mag_phase(H);
STRUCTURE.Sensor(i).H(2).mag = mag;
STRUCTURE.Sensor(i).H(2).phase = phase;


OPTIONS{1} = 'yes';
OPTIONS{2} = 3;
[freqout,H]=tfest(x,f,t,[],OPTIONS);
STRUCTURE.Sensor(i).H(3).total = H;
[mag, phase ] = mag_phase(H);
STRUCTURE.Sensor(i).H(3).mag = mag;
STRUCTURE.Sensor(i).H(3).phase = phase;

STRUCTURE.Sensor(i).Coh = STRUCTURE.Sensor(i).H(1).total./ STRUCTURE.Sensor(i).H(2).total;

STRUCTURE.Sensor(i).FFT = X;


end



STRUCTURE.f = freqout;

end


function [mag, phase ] = mag_phase(H)
mag = 20 * log10(abs(H));
phase = unwrap(angle(H)) * 180 / pi;

end