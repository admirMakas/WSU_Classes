function [ samples ] = SensitivityStudy( X, Number, percent)
%MCS Takes a structured array of inputs with parameters
%   X.distribution = normal, lognormal
%   X.parameters = normal(mean, std) , lognormal(Norm mean, Norm std)

samples = ones(Number*length(X),length(X));

weights = linspace(1-percent, 1+percent, Number)';

for i = 1:length(X)
    mu = X(i).parameters(1);
    samples(1+(Number*(i-1)):Number*i,i) = weights;
    samples(:,i) =  samples(:,i)*mu;

end

end