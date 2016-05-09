%this code uses rand function to generate NSamples for NVar
%the data generated using rand function is the sampled using
%latin hypercube method. Then the data is transformed into normally
%distributed variables. Note rand function is uniformly distributed.
clear all
clc

NSample = 5;
NVar = 2;
dataN=[];

data = rand(NSample,NVar);

    for i = 1:NVar %loop used to work on rows 1 and 2 defined by NVar
        index = randperm(NSample); %permutes the NSamples
        prob = (index'-data(:,i))/NSample; %assuming this is LHS method
        dataN(:,i) = sqrt(2)*erfinv(2*prob-1); %generates normal dist.
    end;
    
hist(data(:,2))
figure
hist(dataN(:,2))