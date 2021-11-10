clear; clc;
% Paper Title: A Generalized Unscented Transformation for Probability Distributions
%
% This code simulates the first set of results for Case Study 1
% It evaluates the quadratic function   3*x + 2*x^2
% It calls the functions that generates the raw moments and central moments
% for the 10 different probability distributions 

addpath('Distribution Moments');
addpath('Unscented Transforms');

%% preallocate for storage
global quad_mean_GenUT quad_mean_HOSPUT quad_mean_UT quad_mean_Monte...
    quad_cov_GenUT quad_cov_HOSPUT quad_cov_UT quad_cov_Monte

quad_mean_GenUT = zeros(10,1); % GenUT approximation of quadratic mean
quad_mean_HOSPUT = zeros(10,1); % HOSPUT approximation of quadratic mean
quad_mean_UT = zeros(10,1); % Scaled UT approximation of quadratic mean
quad_mean_Monte = zeros(10,1); % Monte Carlo approximation of quadratic mean
quad_mean_True = zeros(10,1); % True quadratic mean

quad_cov_GenUT = zeros(10,1); % GenUT approximation of quadratic variance
quad_cov_HOSPUT = zeros(10,1); % HOSPUT approximation of quadratic variance
quad_cov_UT = zeros(10,1); % Scaled UT approximation of quadratic variance
quad_cov_Monte = zeros(10,1); % Monte Carlo approximation of quadratic variance
quad_cov_True = zeros(10,1); % True quadratic variance

%%  Evaluate means and covariances for the 10 different probability distributions
global count_val;
count_val = 0; % Initialize counter as 0
num_Monte = 300000; % Use 300,000 Monte carlo draws

%-------------------------------------------------------------------------%
%********************  Gaussian Random Variable **************************%
%-------------------------------------------------------------------------%

%\***** Gaussian(mu,sigma)         mean = 1,   standard deviation = 2 ********\%
mu = 1;         sigma = 2;
[mu, Gauss_second, Gauss_third, ...
    Gauss_fourth, Gauss_Ex] =  Gaussian_moments(mu, sigma);
% Evalute and store results of the different unscented transforms
Sample_statistics(mu, Gauss_second, Gauss_third, Gauss_fourth ); 
% Evaluate and store Monte Carlo draws
x = normrnd(mu,sigma,1,num_Monte);
y_monte = quadraticTransform(x);
quad_mean_Monte(count_val) = mean(y_monte);
quad_cov_Monte(count_val) = cov(y_monte);
% Evaluate and store the true mean and variance
[quad_mean_True(count_val),quad_cov_True(count_val)] = analyticMean(Gauss_Ex);

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
y_monte = quadraticTransform(x);
quad_mean_Monte(count_val) = mean(y_monte);
quad_cov_Monte(count_val) = cov(y_monte);
% Evaluate and store the true mean and variance
[quad_mean_True(count_val),quad_cov_True(count_val)] = analyticMean(Exp_Ex);


%-------------------------------------------------------------------------%
%**********************  Gamma Random Variable ***************************%
%-------------------------------------------------------------------------%

%\***** Gamma(alpha, beeta)      alpha = 1,  beeta = 2 ********\%
alpha = 1;  % shape parameter
beeta = 2;   % scale parameter
[Gam_mean, Gam_second, Gam_third, ...
    Gam_fourth, Gam_Ex ] =  Gamma_moments(alpha,beeta);
% Evalute and store results of the different unscented transforms
Sample_statistics(Gam_mean, Gam_second, Gam_third, Gam_fourth); 
% Evaluate and store Monte Carlo draws
x = gamrnd(alpha,beeta,1,num_Monte); 
y_monte = quadraticTransform(x);
quad_mean_Monte(count_val) = mean(y_monte);
quad_cov_Monte(count_val) = cov(y_monte);
% Evaluate and store the true mean and variance
[quad_mean_True(count_val),quad_cov_True(count_val)] = analyticMean(Gam_Ex);

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
y_monte = quadraticTransform(x);
quad_mean_Monte(count_val) = mean(y_monte);
quad_cov_Monte(count_val) = cov(y_monte);
% Evaluate and store the true mean and variance
[quad_mean_True(count_val),quad_cov_True(count_val)] = analyticMean(Weib_Ex);

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
y_monte = quadraticTransform(x);
quad_mean_Monte(count_val) = mean(y_monte);
quad_cov_Monte(count_val) = cov(y_monte);
% Evaluate and store the true mean and variance
[quad_mean_True(count_val),quad_cov_True(count_val)] = analyticMean(Rayl_Ex);


%-------------------------------------------------------------------------%
%**********************  beeta Random Variable ***************************%
%-------------------------------------------------------------------------%

%\***** beeta(alpha, beeta)      alpha = 1,  beeta = 2 ********\%
alpha = 1;  % first shape parameter
beeta = 2;   % second shape parameter
[beeta_mean, beeta_second, beeta_third, ...
    beeta_fourth, beeta_Ex ] =  Beta_moments(alpha,beeta);
% Evalute and store results of the different unscented transforms
Sample_statistics(beeta_mean, beeta_second, beeta_third, beeta_fourth); 
% Evaluate and store Monte Carlo draws
x = betarnd(alpha,beeta,1,num_Monte); 
y_monte = quadraticTransform(x);
quad_mean_Monte(count_val) = mean(y_monte);
quad_cov_Monte(count_val) = cov(y_monte);
% Evaluate and store the true mean and variance
[quad_mean_True(count_val),quad_cov_True(count_val)] = analyticMean(beeta_Ex);


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
y_monte = quadraticTransform(x);
quad_mean_Monte(count_val) = mean(y_monte);
quad_cov_Monte(count_val) = cov(y_monte);
% Evaluate and store the true mean and variance
[quad_mean_True(count_val),quad_cov_True(count_val)] = analyticMean(Bino_Ex);

%-------------------------------------------------------------------------%
%*********************  Poisson Random Variable *************************%
%-------------------------------------------------------------------------%

%\***** Poisson(lambda)         lambda = 2 ********\%
lambda = 2;      % rate parameter
[Poiss_mean, Poiss_second, Poiss_third, Poiss_fourth, Poiss_Ex ] =  Poisson_moments(lambda);
% Evalute and store results of the different unscented transforms
Sample_statistics(Poiss_mean, Poiss_second, Poiss_third, Poiss_fourth); 
% Evaluate and store Monte Carlo draws
x = poissrnd(lambda*ones(1, num_Monte)); 
y_monte = quadraticTransform(x);
quad_mean_Monte(count_val) = mean(y_monte);
quad_cov_Monte(count_val) = cov(y_monte);
% Evaluate and store the true mean and variance
[quad_mean_True(count_val),quad_cov_True(count_val)] = analyticMean(Poiss_Ex);

%-------------------------------------------------------------------------%
%********************  Geometric Random Variable *************************%
%-------------------------------------------------------------------------%

%\***** Geometric(p)         p = 0.5 ********\%
p = 0.5;      % probability parameter
[Geo_mean, Geo_second, Geo_third, Geo_fourth, Geo_Ex ] =  Geometric_moments(p);
% Evalute and store results of the different unscented transforms
Sample_statistics(Geo_mean, Geo_second, Geo_third, Geo_fourth); 
% Evaluate and store Monte Carlo draws
x = geornd(p*ones(1, num_Monte)); 
y_monte = quadraticTransform(x);
quad_mean_Monte(count_val) = mean(y_monte);
quad_cov_Monte(count_val) = cov(y_monte);
% Evaluate and store the true mean and variance
[quad_mean_True(count_val),quad_cov_True(count_val)] = analyticMean(Geo_Ex);

%-------------------------------------------------------------------------%
%****************  Negative Binomial Random Variable *********************%
%-------------------------------------------------------------------------%

%\***** NegativeBinomial(r, p)        r = 4,  p = 0.67 ********\%
r = 4;          % number of successes
p = 0.67;       % probability of success
[NB_mean, NB_second, NB_third, NB_fourth, NB_Ex ] =  Negative_Binomial_moments(r,p);
% Evalute and store results of the different unscented transforms
Sample_statistics(NB_mean, NB_second, NB_third, NB_fourth); 
% Evaluate and store Monte Carlo draws
x = nbinrnd(r,p,1, num_Monte); 
y_monte = quadraticTransform(x);
quad_mean_Monte(count_val) = mean(y_monte);
quad_cov_Monte(count_val) = cov(y_monte);
% Evaluate and store the true mean and variance
[quad_mean_True(count_val),quad_cov_True(count_val)] = analyticMean(NB_Ex);

%*************************************************************************%
%********************  Display results as a table ************************%
%*************************************************************************%

%\- For the mean results
% Calculate percentage error
GenUT = abs((quad_mean_GenUT - quad_mean_True)./quad_mean_True)*100;
UT = abs((quad_mean_UT - quad_mean_True)./quad_mean_True)*100;
MC = abs((quad_mean_Monte - quad_mean_True)./quad_mean_True)*100;
HOSPUT = abs((quad_mean_HOSPUT - quad_mean_True)./quad_mean_True)*100;

GenUT = round(GenUT,3);     % Round to 3 decimal places
UT = round(UT,3);           % Round to 3 decimal places
HOSPUT = round(HOSPUT,3);   % Round to 3 decimal places
MC = round(MC,3);   % Round to 3 decimal places

x = {'Gaussian(1,4)'; 'Exponenial(2)'; 'Gamma(1,2)'; 'Weibull(1,2)';...
    'Rayleigh(1)'; 'beeta(3,4)'; 'Binomial(3,0.3)'; 'Poisson(2)';...
    'Geometric(0.5)'; 'NegativeBinomial(4,0.67)'};
Tvar_mean = table(x, GenUT, UT, MC, HOSPUT);
disp('-------------------------------------------------------------------------');
disp('### Percentage error in propagating the mean of  y = 3x + 2x^2 ###');
disp(Tvar_mean);
disp('-------------------------------------------------------------------------');


%\- For the variance results
% Calculate percentage error
GenUT = abs((quad_cov_GenUT - quad_cov_True)./quad_cov_True)*100;
UT = abs((quad_cov_UT - quad_cov_True)./quad_cov_True)*100;
MC = abs((quad_cov_Monte - quad_cov_True)./quad_cov_True)*100;
HOSPUT = abs((quad_cov_HOSPUT - quad_cov_True)./quad_cov_True)*100;

GenUT = round(GenUT,3);     % Round to 3 decimal places
UT = round(UT,3);           % Round to 3 decimal places
HOSPUT = round(HOSPUT,3);   % Round to 3 decimal places
MC = round(MC,3);   % Round to 3 decimal places

x = {'Gaussian(1,4)'; 'Exponenial(2)'; 'Gamma(1,2)'; 'Weibull(1,2)';...
    'Rayleigh(1)'; 'beeta(3,4)'; 'Binomial(3,0.3)'; 'Poisson(2)';...
    'Geometric(0.5)'; 'NegativeBinomial(4,0.67)'};
Tvar_cov = table(x, GenUT, UT, MC, HOSPUT);
fprintf(' \n ');
fprintf(' \n ');
disp('-------------------------------------------------------------------------');
disp('### Percentage error in propagating the variance of  y = 3x + 2x^2 ###');
disp(Tvar_cov);
disp('-------------------------------------------------------------------------');



%%  Nonlinear quadratic function

function y = quadraticTransform(x)
y = 3*x + 2*x.^2;
end

%% Function that calculates the analytic mean and variance of 3x + 2x^2
function [y_mean_true,y_cov_true] = analyticMean(Ex)
Ex1 = Ex(1);    % E [ x ]  
Ex2 = Ex(2);    % E [ x^2 ]
Ex3 = Ex(3);    % E [ x^3 ]
Ex4 = Ex(4);    % E [ x^4 ]

y_mean_true = 3*Ex1 + 2*Ex2;
y_cov_true = 4*Ex4 + 12*Ex3 +9*Ex2 -2*y_mean_true*(3*Ex1 + 2*Ex2) + y_mean_true^2;
end

%% Function that evaluates and stores the approximate means and covariances
function Sample_statistics(x_mu, x_cov, x_skew, x_kurt )
% INPUTS
%
%   x_mu        -   Mean of random variable
%   x_cov       -   Variance of random variable
%   x_skew      -   Skewness of random variable
%   x_kurt      -   Kurtosis of random variable

global quad_mean_GenUT quad_mean_HOSPUT quad_mean_UT ...
    quad_cov_GenUT quad_cov_HOSPUT quad_cov_UT count_val
count_val = count_val+1;
% Generate the ensembles and weights for the different transforms
[sigma_UT, weights_UT] = unscentedEnsemble(x_mu, x_cov, sqrt(3)); % scaled UT
[sigma_GenUT, weights_GenUT] = GenUT_Ensemble(x_mu, x_cov, x_skew, x_kurt); % GenUT
[sigma_HOSPUT, weights_HOSPUT] = HOSP(x_mu, x_cov, x_skew, x_kurt); % HOSPUT

% Evaluate and store scaled UT mean and variance
[quad_mean_UT(count_val),quad_cov_UT(count_val),~, ~] = Evaluate_sample_statistics(quadraticTransform(sigma_UT),weights_UT);

% Evaluate and store GenUT mean and variance
[quad_mean_GenUT(count_val),quad_cov_GenUT(count_val),~, ~] = Evaluate_sample_statistics(quadraticTransform(sigma_GenUT),weights_GenUT);

% Evaluate and store HOSPUT mean and variance
[quad_mean_HOSPUT(count_val),quad_cov_HOSPUT(count_val),~, ~] = Evaluate_sample_statistics(quadraticTransform(sigma_HOSPUT),weights_HOSPUT);
end

