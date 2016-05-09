function [ g, gDelta ] = LSB( X )
% Failure event is defined as $delta_{max} > 2.0$
Failure = 2;

g = Failure - delta_max(X(:,1), X(:,2), X(:,3));  % this is the LSB
gDelta = delta_delta_max(X(:,1), X(:,2), X(:,3)); % Gradients of the LSB

end

function [ delta ] = delta_max(P,E,W)
L = 30*12;                  % was ft now in
I = 1.33*10^3;              % in^4

delta =  (P.*L.^3)./ (48.*E.*I) + (5*W.*L.^4)./(385.*E.*I);

end

function [ delta_delta ] = delta_delta_max(P,E,W)
L = 30*12;                  % was ft now in
I = 1.33*10^3;              % in^4

dmax_dP = -(L.^3)./ (48.*E.*I);
dmax_dE = (P.*L.^3)./ (48.*E.^2.*I) + (5*W.*L.^4)./(385.*E.^2.*I);
dmax_dW = -(5*L.^4)./(385.*E.*I);

delta_delta = [dmax_dP; dmax_dE; dmax_dW];

end