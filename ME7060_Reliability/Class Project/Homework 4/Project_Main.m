%% Project                                       Student:      Daniel Clark
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

MCS_samples = MCS(X, 10 );
LHS_samples = MCS_LHS(X, 10);
Sensitivity_samples = SensitivityStudy(X, 10, 0.3);

%% Adding folder to Path 
CasePath = 'C:\Temp'; addpath(CasePath,'-end');
specificRun = Sensitivity_samples;
[runs,~] = size(specificRun);

%% Running Analysis
freqArray = zeros(runs,3);



tic
for i = 1:1
    
    % Change Information for this run
    Change2 = specificRun(i,:);
    
    newFileName = strcat('Job_temp','.inp');
    newFilePath = sprintf('%s\\%s',CasePath,newFileName);
    GenerateNewJob('Job_original.dat', newFilePath, Change2, newFileName)
    
    
%     JOB = 'Job_temp';
%     Enter =  sprintf('abaqus job=%s inp=%s interactive', JOB, newFileName);
%     dos(Enter)
%     
%     frequencies = ReadFile_dat( 'Job_temp.dat' );
%     displacementArray(i,:) = frequencies;
%     save('solution.mat','PlaceHolders', 'displacementArray');
end
toc


% Average = mean(abs(displacementArray))
% Stdev = std(abs(displacementArray))
% nf = sum(abs(displacementArray)>2.0)
% Pf = nf/length(displacementArray)



