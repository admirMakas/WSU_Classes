function [ MCS_samples ] = MCS( X, Number_of_runs )
%MCS Takes a structured array of inputs with parameters
%   X.distribution = normal, lognormal
%   X.parameters = normal(mean, std) , lognormal(Norm mean, Norm std)

MCS_samples = zeros(Number_of_runs,length(X));

for i = 1:length(X)
    
    if sum(X(i).distribution == 'normal000') == 9
        mu = X(i).parameters(1);
        sigma = X(i).parameters(2);
        MCS_samples(:,i) = mu + sigma.*icdf('normal',...
            rand(Number_of_runs,1),0,1);
    
    elseif sum(X(i).distribution == 'lognormal') == 9
        mu = X(i).parameters(1);
        sigma = X(i).parameters(2);
        mu_log = log((mu^2)/sqrt(sigma^2+mu^2));
        sigma_log = sqrt(log((sigma^2/(mu^2))+1));
        MCS_samples(:,i) = exp(mu_log + sigma_log.*icdf('normal',...
            rand(Number_of_runs,1),0,1));
        
    else
        T = sprintf('Not a defined distribution type');
        error(T)
        
    end
    
end

end