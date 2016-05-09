%% Sensitivity Study                             Student:      Daniel Clark
clear all                                      % Instructor: Dr. Ha-Rok Bae
close all                                      % Class: ME 7060 Spring 2016
clc
format shorte, warning('off')

%% My system
mu_Kz = 50*10^6;             
sigma_Kz = 0;

mu_T = 1;
sigma_T = 0;

mu_L = 1;
sigma_L = 0;

mu_P = 0.33;
sigma_P = 0;

mu_E = 68.98 * 10^6;
sigma_E = 0;

%% Distributions and Parameters
X(1).distribution = 'normal000'; X(1).parameters = [mu_Kz, sigma_Kz]; 
X(2).distribution = 'normal000'; X(2).parameters = [mu_T, sigma_T];
X(3).distribution = 'normal000'; X(3).parameters = [mu_L, sigma_L];
X(4).distribution = 'normal000'; X(4).parameters = [mu_P, sigma_P];
X(5).distribution = 'lognormal'; X(5).parameters = [mu_E, sigma_E];

%% Samples
Number = 10;
percent = 0.3;
Sensitivity_samples = SensitivityStudy( X, Number, percent);

%% Adding folder to Path 
CasePath = 'C:\Temp'; 
nastranPath = 'C:\MSC.Software\MSC_Nastran\20130\bin\nastran.exe';
%CasePath = 'C:\Users\ecslogon\Google Drive\Reliability\Project';
addpath(CasePath,'-end');



specificRun = Sensitivity_samples;
[runs,~] = size(specificRun);

%% Running Analysis
freqArray = zeros(runs,3);

tic
for i = 1    
    % Change Information for this run
    Change2 = specificRun(i,:);
    
    newFileName = strcat('Job_temp','.dat');
    newFilePath = sprintf('%s\\%s',CasePath,newFileName);
    GenerateNewJob('Job_original.dat', newFilePath, Change2, newFileName)
    
    JOB = 'Job_temp.dat';
    
    %Enter =  sprintf('%s\\%s', CasePath, JOB);
    
    Enter =  sprintf('nastran %s bat=no mem=4 parallel=2' , JOB);
    
    %dos(nastranPath)
    dos(Enter)
    
    frequencies = ReadFile( 'Job_temp.f06' );
    freqArray(i,:) = frequencies;
    %save('solution.mat','specificRun', 'freqArray');
end
toc
break

for i = 1:length(X)
    figure
    X = Sensitivity_samples(1+(Number*(i-1)):Number*i,i);
    Y = freqArray(1+(Number*(i-1)):Number*i,1);
    plot(X,Y)
end





