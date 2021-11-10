% Paper Title: A Generalized Unscented Transformation for Probability Distributions
%
function [mu, second_cen_moment, third_cen_moment, ...
    fourth_cen_moment, Ex] =  Geometric_moments(p)
% INPUTS
%
%   p                       -   rate parameter
%
% OUTPUTS
%
%   mu                      -   Mean of random variable
%   second_cen_moment       -   Variance of random variable
%   third_cen_moment        -   Skewness of random variable
%   fourth_cen_moment       -   Kurtosis of random variable
%   Ex                      -   Vector of first four raw mments


% syms p t
syms t
J = (p /  (1 - (1-p)*exp(t) )  ); % MGF
M1 = diff(J,t);
M2 = diff(M1,t);
M3 = diff(M2,t);
M4 = diff(M3,t);

Ex1 = subs(M1, t, 0);     % E [ x ]
Ex2 = subs(M2, t, 0);     % E [ x^2 ]
Ex3 = subs(M3, t, 0);     % E [ x^3 ]
Ex4 = subs(M4, t, 0);     % E [ x^4 ]

% The mean
mu = Ex1;    % = ( (1 - p))/p
% Variance  --->   E[ (x - mu)^2 ]
second_cen_moment = simplify(Ex2 - mu*mu);  % ( (1 - p))/p^2
% Skew      --->   E[ (x - mu)^3 ]
third_cen_moment = simplify(Ex3 - 3*Ex2*mu + 3*mu*mu^2 - mu^3); % (p^2 - 3*p + 2)/p^3
% Kurtosis  --->   E[ (x - mu)^4 ]
fourth_cen_moment = simplify(Ex4 - 4*Ex3*mu + 6*Ex2*mu^2 - 4*mu*mu^3 + mu^4); % (1-p)*( p^2  - 9*p  + 9)/p^4

%% Save output values as double
Ex1 = double(Ex1);      Ex2 = double(Ex2);
Ex3 = double(Ex3);      Ex4 = double(Ex4);

Ex = [Ex1; Ex2; Ex3; Ex4];
mu = double(mu);             
second_cen_moment = double(second_cen_moment);
third_cen_moment = double(third_cen_moment);
fourth_cen_moment = double(fourth_cen_moment);

end