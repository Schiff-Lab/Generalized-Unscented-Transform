% Paper Title: A Generalized Unscented Transformation for Probability Distributions
%
function [mu, second_cen_moment, third_cen_moment, ...
    fourth_cen_moment, Ex] =  Binomial_moments(n, p)
% INPUTS
%
%   n                       -   number of trials
%   p                       -   probability of success
%
% OUTPUTS
%
%   mu                      -   Mean of random variable
%   second_cen_moment       -   Variance of random variable
%   third_cen_moment        -   Skewness of random variable
%   fourth_cen_moment       -   Kurtosis of random variable
%   Ex                      -   Vector of first four raw mments


% syms n p t
syms t
J = ( 1-p + p*exp(t))^n; % MGF
M1 = diff(J,t);
M2 = diff(M1,t);
M3 = diff(M2,t);
M4 = diff(M3,t);

Ex1 = subs(M1, t, 0);     % E [ x ]
Ex2 = subs(M2, t, 0);     % E [ x^2 ]
Ex3 = subs(M3, t, 0);     % E [ x^3 ]
Ex4 = subs(M4, t, 0);     % E [ x^4 ]

% The mean
mu = Ex1;    % = np
% Variance  --->   E[ (x - mu)^2 ]
second_cen_moment = simplify(Ex2 - mu*mu);  % np(1-p)
% Skew      --->   E[ (x - mu)^3 ]
third_cen_moment = simplify(Ex3 - 3*Ex2*mu + 3*mu*mu^2 - mu^3); % n*p*(2*p^2 - 3*p + 1)
% Kurtosis  --->   E[ (x - mu)^4 ]
fourth_cen_moment = simplify(Ex4 - 4*Ex3*mu + 6*Ex2*mu^2 - 4*mu*mu^3 + mu^4); % n*p*(1-p)*(3*n*p - 6*p - 3*n*p^2 + 6*p^2 + 1)

%% Save output values as double
Ex1 = double(Ex1);      Ex2 = double(Ex2);
Ex3 = double(Ex3);      Ex4 = double(Ex4);

Ex = [Ex1; Ex2; Ex3; Ex4];
mu = double(mu);             
second_cen_moment = double(second_cen_moment);
third_cen_moment = double(third_cen_moment);
fourth_cen_moment = double(fourth_cen_moment);

end