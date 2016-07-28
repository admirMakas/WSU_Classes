close all
clear all
clc

load('MatrixDeflation.mat')
start = 29; % Removing Z rx ry

while start > 1
    ints = start-2:start;
    start = start-6;
    Kfull(:,ints) = [];
    Kfull(ints,:) = [];
    Mfull(:,ints) = [];
    Mfull(ints,:) = [];
end

BC = [1,2]; REDUCED_K = Kfull; REDUCED_M = Mfull;

for i = linspace(2,1,2)
   REDUCED_K(:,BC(i)) = []; % col
   REDUCED_M(:,BC(i)) = []; % col
   REDUCED_K(BC(i),:) = []; % row
   REDUCED_M(BC(i),:) = []; % row
end

%A = REDUCED_M\REDUCED_K;
clear all
clc

A = [1,2,3;2,5,7;3,7,8]; % 

[ModeShapes,Freq_Sqrt] = eig(A);

Freq_Sqrt = diag(Freq_Sqrt);
[~,idx] = sort(Freq_Sqrt);
Freq_Sqrt = Freq_Sqrt(idx);
ModeShapes = ModeShapes(:,idx);
num = 2;
B = A - Freq_Sqrt(num)*ModeShapes(:,num)*ModeShapes(:,num)';

[ModeShapes1,Freq_Sqrt1] = eig(B);




break
thing = zeros(3,3);

for i = 2:3
    thing = thing+ Freq_Sqrt(i)*ModeShapes(:,i)*ModeShapes(:,i)';
end