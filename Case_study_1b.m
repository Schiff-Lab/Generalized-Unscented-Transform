clear; clc;
% Paper Title: A Generalized Unscented Transformation for Probability Distributions
%
% This code simulates the second set of results for Case Study 1
% It evaluates the trigonometric function   y = sin(x)
% It calls the functions that generates the raw moments and central moments
% for the 10 different probability distributions 

addpath('Distribution Moments');
addpath('Unscented Transforms');
%% preallocate for storage
global trig_mean_GenUT trig_mean_HOSPUT trig_mean_UT trig_mean_Monte...
    trig_cov_GenUT trig_cov_HOSPUT trig_cov_UT trig_cov_Monte

trig_mean_GenUT = zeros(10,1); % GenUT approximation of trigonometric mean
trig_mean_HOSPUT = zeros(10,1); % HOSPUT approximation of trigonometric mean
trig_mean_UT = zeros(10,1); % Scaled UT approximation of trigonometric mean
trig_mean_Monte = zeros(10,1); % Monte Carlo approximation of trigonometric mean
trig_mean_True = zeros(10,1); % True quadratic mean

trig_cov_GenUT = zeros(10,1); % GenUT approximation of trigonometric variance
trig_cov_HOSPUT = zeros(10,1); % HOSPUT approximation of trigonometric variance
trig_cov_UT = zeros(10,1); % Scaled UT approximation of trigonometric variance
trig_cov_Monte = zeros(10,1); % Monte Carlo approximation of trigonometric variance
trig_cov_True = zeros(10,1); % True trigonometric variance

%%  Evaluate means and covariances for the 10 different probability distributions
global count_val;
count_val = 0; % Initialize counter as 0
num_Monte = 300000; % Use 300,000 Monte carlo draws

%-------------------------------------------------------------------------%
%********************  Gaussian Random Variable **************************%
%-------------------------------------------------------------------------%

%\*** Gaussian(mu,sigma)    mean = 1.57,  standard deviation = 0.3162 ***\%
mu = 1.57;         sigma = sqrt(0.1);
[mu, Gauss_second, Gauss_third, ...
    Gauss_fourth, Gauss_Ex] =  Gaussian_moments(mu, sigma);
% Evalute and store results of the different unscented transforms
Sample_statistics(mu, Gauss_second, Gauss_third, Gauss_fourth ); 
% Evaluate and store Monte Carlo draws
x = normrnd(mu,sigma,1,num_Monte);
y_monte = trigonometricTransform(x);
trig_mean_Monte(count_val) = mean(y_monte);
trig_cov_Monte(count_val) = cov(y_monte);
% Evaluate and store the true mean and variance
E_sin_theta = exp(-sigma^2/2)*sin(mu);       % E [  sin(theta) ]
E_cos_2_theta = exp(-2*sigma^2)*cos(2*mu);   % E [  cos(2*theta) ]
[trig_mean_True(count_val),trig_cov_True(count_val)] = analyticMean(E_sin_theta, E_cos_2_theta);

%-------------------------------------------------------------------------%
%*******************  Exponential Random Variable ************************%
%-------------------------------------------------------------------------%

%\***** Exponential(lambda)      lambda = 2 ********\%
lambda = 2;
[Exp_mean, Exp_second, Exp_third, ...
    Exp_fourth, Exp_Ex ] =  Exponential_moments(lambda);
% Evalute and store results of the different unscented transforms
Sample_statistics(Exp_mean, Exp_second, Exp_third, Exp_fourth); 
% Evaluate and store Monte Carlo draws
x = exprnd((1/lambda),1,num_Monte); % Adjust to account for Matlab's way of its PMF representation
                                 % 1/lambda = mu. Matlab uses mu
y_monte = trigonometricTransform(x);
trig_mean_Monte(count_val) = mean(y_monte);
trig_cov_Monte(count_val) = cov(y_monte);
% Evaluate and store the true mean and variance
E_sin_theta = (lambda)/(lambda^2 + 1);  % E [  sin(theta) ]
E_cos_2_theta = (lambda^2)/(lambda^2 + 4);  % E [  cos(2*theta) ]
[trig_mean_True(count_val),trig_cov_True(count_val)] = analyticMean(E_sin_theta, E_cos_2_theta);


%-------------------------------------------------------------------------%
%**********************  Gamma Random Variable ***************************%
%-------------------------------------------------------------------------%

%\***** Gamma(alpha, beeta)      alpha = 0.5,  beeta = 0.5 ********\%
alpha = 0.5;  % shape parameter
beeta = 0.5;   % scale parameter
[Gam_mean, Gam_second, Gam_third, ...
    Gam_fourth, Gam_Ex ] =  Gamma_moments(alpha,beeta);
% Evalute and store results of the different unscented transforms
Sample_statistics(Gam_mean, Gam_second, Gam_third, Gam_fourth); 
% Evaluate and store Monte Carlo draws
x = gamrnd(alpha,beeta,1,num_Monte); 
y_monte = trigonometricTransform(x);
trig_mean_Monte(count_val) = mean(y_monte);
trig_cov_Monte(count_val) = cov(y_monte);
% Evaluate and store the true mean and variance
E_sin_theta = sin(alpha*atan(beeta)) / ( sqrt(1+ beeta^2) )^alpha; % E [  sin(theta) ]
E_cos_2_theta = cos(-alpha*atan(2*beeta)) / ( sqrt(1+ 2^2*beeta^2) )^alpha;  % E [  cos(2*theta) ]
[trig_mean_True(count_val),trig_cov_True(count_val)] = analyticMean(E_sin_theta, E_cos_2_theta);

%-------------------------------------------------------------------------%
%*********************  Weibull Random Variable **************************%
%-------------------------------------------------------------------------%

%\***** Weibull(alpha, beeta)      alpha = 1,  beeta = 2 ********\%
alpha = 1;  % scale parameter
beeta = 2;   % shape parameter

[Weib_mean, Weib_second, Weib_third, ...
    Weib_fourth, Weib_Ex ] =  Weibull_moments(alpha,beeta);
% Evalute and store results of the different unscented transforms
Sample_statistics(Weib_mean, Weib_second, Weib_third, Weib_fourth); 
% Evaluate and store Monte Carlo draws
x = wblrnd(alpha,beeta,1,num_Monte); 
y_monte = trigonometricTransform(x);
trig_mean_Monte(count_val) = mean(y_monte);
trig_cov_Monte(count_val) = cov(y_monte);
% Evaluate and store the true mean and variance
syms k 
t = 1;
E_sin_theta = symsum(gamma(1 + k/beeta)*(1i*t)^k*alpha^k/factorial(k) ,k,0,Inf)     ;
E_sin_theta = imag(  double(E_sin_theta)  ); % E [  sin(theta) ]
t = 2;
E_cos_2_theta = symsum(gamma(1 + k/beeta)*(1i*t)^k*alpha^k/factorial(k) ,k,0,Inf) ;
E_cos_2_theta = real(  double(E_cos_2_theta)  ); % E [  cos(2*theta) ]
[trig_mean_True(count_val),trig_cov_True(count_val)] = analyticMean(E_sin_theta, E_cos_2_theta);
clear k t

%-------------------------------------------------------------------------%
%*********************  Rayleigh Random Variable *************************%
%-------------------------------------------------------------------------%

%\***** Rayleigh(sigma)      sigma = 1   ********\%
sigma = 1;  % scale parameter

[Rayl_mean, Rayl_second, Rayl_third, ...
    Rayl_fourth, Rayl_Ex ] =  Rayleigh_moments(sigma);
% Evalute and store results of the different unscented transforms
Sample_statistics(Rayl_mean, Rayl_second, Rayl_third, Rayl_fourth); 
% Evaluate and store Monte Carlo draws
x = raylrnd(sigma*ones(1,num_Monte)); 
y_monte = trigonometricTransform(x);
trig_mean_Monte(count_val) = mean(y_monte);
trig_cov_Monte(count_val) = cov(y_monte);
% Evaluate and store the true mean and variance
t = 1;
E_sin_theta = 1 - sigma*t*exp(-sigma^2*t^2/2)*sqrt(pi/2)*( erfi(sigma*t/sqrt(2))   - 1i )     ;
E_sin_theta = imag(  double(E_sin_theta)  );    % E [  sin(theta) ]
t = 2;
E_cos_2_theta = 1 - sigma*t*exp(-sigma^2*t^2/2)*sqrt(pi/2)*( erfi(sigma*t/sqrt(2))   - 1i )  ;
E_cos_2_theta = real(  double(E_cos_2_theta)  );    % E [  cos(2*theta) ]
[trig_mean_True(count_val),trig_cov_True(count_val)] = analyticMean(E_sin_theta, E_cos_2_theta);
clear t

%-------------------------------------------------------------------------%
%**********************  Beta Random Variable ***************************%
%-------------------------------------------------------------------------%

%\***** beeta(alpha, beeta)      alpha = 3,  beeta = 4 ********\%
alpha = 3;  % first shape parameter
beeta = 4;   % second shape parameter
[beeta_mean, beeta_second, beeta_third, ...
    beeta_fourth, beeta_Ex ] =  Beta_moments(alpha,beeta);
% Evalute and store results of the different unscented transforms
Sample_statistics(beeta_mean, beeta_second, beeta_third, beeta_fourth); 
% Evaluate and store Monte Carlo draws
x = betarnd(alpha,beeta,1,num_Monte); 
y_monte = trigonometricTransform(x);
trig_mean_Monte(count_val) = mean(y_monte);
trig_cov_Monte(count_val) = cov(y_monte);
% Evaluate and store the true mean and variance
syms k
t = 1;
E_sin_theta = symsum(( beta(alpha+k,beeta)/ beta(alpha,beeta)   )*(1i*t)^k/factorial(k) ,k,0,Inf)     ;
E_sin_theta = imag(  double(E_sin_theta)  );    % E [  sin(theta) ]
t = 2;
E_cos_2_theta = symsum(( beta(alpha+k,beeta)/ beta(alpha,beeta)   )*(1i*t)^k/factorial(k) ,k,0,Inf) ;
E_cos_2_theta = real(  double(E_cos_2_theta)  );    % E [  cos(2*theta) ]
[trig_mean_True(count_val),trig_cov_True(count_val)] = analyticMean(E_sin_theta, E_cos_2_theta);
clear t k

%-------------------------------------------------------------------------%
%*********************  Binomial Random Variable *************************%
%-------------------------------------------------------------------------%

%\***** Binomial(n, p)      n = 3,  p = 0.3 ********\%
n = 3;      % number of trials
p = 0.3;    % probability of success
[Bino_mean, Bino_second, Bino_third, ...
    Bino_fourth, Bino_Ex ] =  Binomial_moments(n,p);
% Evalute and store results of the different unscented transforms
Sample_statistics(Bino_mean, Bino_second, Bino_third, Bino_fourth); 
% Evaluate and store Monte Carlo draws
x = binornd(n,p,1,num_Monte); 
y_monte = trigonometricTransform(x);
trig_mean_Monte(count_val) = mean(y_monte);
trig_cov_Monte(count_val) = cov(y_monte);
% Evaluate and store the true mean and variance
x = 1 - p;   
E_sin_theta = 0;  E_cos_2_theta = 0;
for k = 0:n
    y = cos(k*2) ; 
    E_cos_2_theta = E_cos_2_theta + nchoosek(n,k)*p^k*y*x^(n-k); % E [  cos(2*theta) ]
    y =  sin(k*1); 
    E_sin_theta = E_sin_theta + nchoosek(n,k)*p^k*y*x^(n-k);    % E [  sin(theta) ]
end
[trig_mean_True(count_val),trig_cov_True(count_val)] = analyticMean(E_sin_theta, E_cos_2_theta);


%-------------------------------------------------------------------------%
%*********************  Poisson Random Variable *************************%
%-------------------------------------------------------------------------%

%\***** Poisson(lambda)         lambda = 0.1 ********\%
lambda = 0.1;      % rate parameter
[Poiss_mean, Poiss_second, Poiss_third, Poiss_fourth, Poiss_Ex ] =  Poisson_moments(lambda);
% Evalute and store results of the different unscented transforms
Sample_statistics(Poiss_mean, Poiss_second, Poiss_third, Poiss_fourth); 
% Evaluate and store Monte Carlo draws
x = poissrnd(lambda*ones(1, num_Monte)); 
y_monte = trigonometricTransform(x);
trig_mean_Monte(count_val) = mean(y_monte);
trig_cov_Monte(count_val) = cov(y_monte);
% Evaluate and store the true mean and variance
E_sin_theta = exp(lambda*(cos(1) - 1))*sin(lambda*sin(1));   % E [  sin(theta) ]
E_cos_2_theta = exp(lambda*(cos(2) - 1))*cos(lambda*sin(2)); % E [  cos(2*theta) ]
[trig_mean_True(count_val),trig_cov_True(count_val)] = analyticMean(E_sin_theta, E_cos_2_theta);

%-------------------------------------------------------------------------%
%********************  Geometric Random Variable *************************%
%-------------------------------------------------------------------------%

%\***** Geometric(p)         p = 0.7 ********\%
p = 0.7;      % probability parameter
[Geo_mean, Geo_second, Geo_third, Geo_fourth, Geo_Ex ] =  Geometric_moments(p);
% Evalute and store results of the different unscented transforms
Sample_statistics(Geo_mean, Geo_second, Geo_third, Geo_fourth); 
% Evaluate and store Monte Carlo draws
x = geornd(p*ones(1, num_Monte)); 
y_monte = trigonometricTransform(x);
trig_mean_Monte(count_val) = mean(y_monte);
trig_cov_Monte(count_val) = cov(y_monte);
% Evaluate and store the true mean and variance
q = 1-p;
E_sin_theta =  p*q*sin(1)/(  (1 - q*cos(1))^2 + q^2*(sin(1))^2 );   % E [  sin(theta) ]
E_cos_2_theta = p*(1-q*cos(2))/(  (1 - q*cos(2))^2 + q^2*(sin(2))^2 );  % E [  cos(2*theta) ]
[trig_mean_True(count_val),trig_cov_True(count_val)] = analyticMean(E_sin_theta, E_cos_2_theta);

%-------------------------------------------------------------------------%
%****************  Negative Binomial Random Variable *********************%
%-------------------------------------------------------------------------%

%\***** NegativeBinomial(r, p)        r = 4,  p = 0.67 ********\%
r = 0.4;          % number of successes
p = 0.67;       % probability of success
[NB_mean, NB_second, NB_third, NB_fourth, NB_Ex ] =  Negative_Binomial_moments(r,p);
% Evalute and store results of the different unscented transforms
Sample_statistics(NB_mean, NB_second, NB_third, NB_fourth); 
% Evaluate and store Monte Carlo draws
x = nbinrnd(r,p,1, num_Monte); 
y_monte = trigonometricTransform(x);
trig_mean_Monte(count_val) = mean(y_monte);
trig_cov_Monte(count_val) = cov(y_monte);
% Evaluate and store the true mean and variance
E_sin_theta =  imag( double(  (p /  (1 - (1-p)*exp(1i*1) )  )^r      ) ); % E [  sin(theta) ]
E_cos_2_theta = real( double(    (p /  (1 - (1-p)*exp(1i*2) )  )^r  ) );  % % E [  cos(2*theta) ]
[trig_mean_True(count_val),trig_cov_True(count_val)] = analyticMean(E_sin_theta, E_cos_2_theta);

%*************************************************************************%
%********************  Display results as a table ************************%
%*************************************************************************%

%\- For the mean results
% Calculate percentage error
GenUT = abs((trig_mean_GenUT - trig_mean_True)./trig_mean_True)*100;
UT = abs((trig_mean_UT - trig_mean_True)./trig_mean_True)*100;
MC = abs((trig_mean_Monte - trig_mean_True)./trig_mean_True)*100;
HOSPUT = abs((trig_mean_HOSPUT - trig_mean_True)./trig_mean_True)*100;

GenUT = round(GenUT,3);     % Round to 3 decimal places
UT = round(UT,3);           % Round to 3 decimal places
HOSPUT = round(HOSPUT,3);   % Round to 3 decimal places
MC = round(MC,3);   % Round to 3 decimal places

x = {'Gaussian(1,4)'; 'Exponenial(2)'; 'Gamma(1,2)'; 'Weibull(1,2)';...
    'Rayleigh(1)'; 'beeta(3,4)'; 'Binomial(3,0.3)'; 'Poisson(2)';...
    'Geometric(0.5)'; 'NegativeBinomial(0.4,0.67)'};
Tvar_mean = table(x, GenUT, UT, MC, HOSPUT);
disp('-------------------------------------------------------------------------');
disp('### Percentage error in propagating the mean of  y = sin(x) ###');
disp(Tvar_mean);
disp('-------------------------------------------------------------------------');


%\- For the variance results
% Calculate percentage error
GenUT = abs((trig_cov_GenUT - trig_cov_True)./trig_cov_True)*100;
UT = abs((trig_cov_UT - trig_cov_True)./trig_cov_True)*100;
MC = abs((trig_cov_Monte - trig_cov_True)./trig_cov_True)*100;
HOSPUT = abs((trig_cov_HOSPUT - trig_cov_True)./trig_cov_True)*100;

GenUT = round(GenUT,3);     % Round to 3 decimal places
UT = round(UT,3);           % Round to 3 decimal places
HOSPUT = round(HOSPUT,3);   % Round to 3 decimal places
MC = round(MC,3);   % Round to 3 decimal places

x = {'Gaussian(1.57,0.1)'; 'Exponenial(2)'; 'Gamma(0.5,0.5)'; 'Weibull(1,2)';...
    'Rayleigh(1)'; 'beeta(3,4)'; 'Binomial(3,0.3)'; 'Poisson(0.1)';...
    'Geometric(0.7)'; 'NegativeBinomial(4,0.67)'};
Tvar_cov = table(x, GenUT, UT, MC, HOSPUT);
fprintf(' \n ');
fprintf(' \n ');
disp('-------------------------------------------------------------------------');
disp('### Percentage error in propagating the variance of  y = sin(x) ###');
disp(Tvar_cov);
disp('-------------------------------------------------------------------------');



%%  Nonlinear quadratic function

function y = trigonometricTransform(x)
y = sin(x);
end

%% Function that calculates the analytic mean and variance of sin(x)
function [y_mean_true,y_cov_true] = analyticMean(E_sin_theta,E_cos_2_theta)
%
%   INPUTS
%
%   E_sin_theta     -   E [ sin(x) ] 
%   E_cos_2_theta   -   E [ cos(2x) ]

y_mean_true = E_sin_theta;
y_cov_true = 0.5*(1 - E_cos_2_theta) - y_mean_true^2; % E[ (sin(x) - E [sin(x)] )^2]
end

%% Function that evaluates and stores the approximate means and covariances
function Sample_statistics(x_mu, x_cov, x_skew, x_kurt )
% INPUTS
%
%   x_mu        -   Mean of random variable
%   x_cov       -   Variance of random variable
%   x_skew      -   Skewness of random variable
%   x_kurt      -   Kurtosis of random variable

global trig_mean_GenUT trig_mean_HOSPUT trig_mean_UT ...
    trig_cov_GenUT trig_cov_HOSPUT trig_cov_UT count_val
count_val = count_val+1;
% Generate the ensembles and weights for the different transforms
[sigma_UT, weights_UT] = unscentedEnsemble(x_mu, x_cov, sqrt(3)); % scaled UT
[sigma_GenUT, weights_GenUT] = GenUT_Ensemble(x_mu, x_cov, x_skew, x_kurt); % GenUT
[sigma_HOSPUT, weights_HOSPUT] = HOSP(x_mu, x_cov, x_skew, x_kurt); % HOSPUT

% Evaluate and store scaled UT mean and variance
[trig_mean_UT(count_val),trig_cov_UT(count_val),~, ~] = Evaluate_sample_statistics(trigonometricTransform(sigma_UT),weights_UT);

% Evaluate and store GenUT mean and variance
[trig_mean_GenUT(count_val),trig_cov_GenUT(count_val),~, ~] = Evaluate_sample_statistics(trigonometricTransform(sigma_GenUT),weights_GenUT);

% Evaluate and store HOSPUT mean and variance
[trig_mean_HOSPUT(count_val),trig_cov_HOSPUT(count_val),~, ~] = Evaluate_sample_statistics(trigonometricTransform(sigma_HOSPUT),weights_HOSPUT);
end

