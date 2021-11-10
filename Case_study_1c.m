clear; clc;
% Paper Title: A Generalized Unscented Transformation for Probability Distributions
%
% This code simulates the third set of results for Case Study 1
% It evaluates the trigonometric function   y = sin(x)
% It calls the functions that generates the raw moments and central moments
% for the Weibull and Poisson distributions 

%% preallocate for storage
global trig_mean_GenUT  trig_mean_UT ...
    trig_cov_GenUT  trig_cov_UT 

%%  Evaluate means and covariances for Weibull and Poisson distributions


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
% Evaluate and store the true mean and variance
syms k 
t = 1;
E_sin_theta = symsum(gamma(1 + k/beeta)*(1i*t)^k*alpha^k/factorial(k) ,k,0,Inf)     ;
E_sin_theta = imag(  double(E_sin_theta)  ); % E [  sin(theta) ]
t = 2;
E_cos_2_theta = symsum(gamma(1 + k/beeta)*(1i*t)^k*alpha^k/factorial(k) ,k,0,Inf) ;
E_cos_2_theta = real(  double(E_cos_2_theta)  ); % E [  cos(2*theta) ]
[trig_mean_True ,trig_cov_True ] = analyticMean(E_sin_theta, E_cos_2_theta);
clear k t

% Evaluate and store Monte Carlo draws
num_iter = 20;
Method = ["300 Monte"; "25000 Monte"; "100000 Monte"; "400000 Monte"; "GenUT" ; "UT"];
Mean_per = zeros(6,num_iter);   % Preallocate for storage
Cov_per = zeros(6,num_iter);    % Preallocate for storage
for i = 1:num_iter
% 300 Monte Carlo Draws
x = wblrnd(alpha,beeta,1, 300);
y_monte300 = trigonometricTransform(x);
y_monte_mean300 = mean(y_monte300);
y_monte_cov300 = cov(y_monte300);

% 25000 Monte Carlo Draws
x = wblrnd(alpha,beeta,1, 25000);
y_monte25000 = trigonometricTransform(x);
y_monte_mean25000 = mean(y_monte25000);
y_monte_cov25000 = cov(y_monte25000);

% 100000 Monte Carlo Draws
x = wblrnd(alpha,beeta,1, 100000);
y_monte100000 = trigonometricTransform(x);
y_monte_mean100000 = mean(y_monte100000);
y_monte_cov100000 = cov(y_monte100000);

% 400000 Monte Carlo Draws
x = wblrnd(alpha,beeta,1, 400000);
y_monte400000 = trigonometricTransform(x);
y_monte_mean400000 = mean(y_monte400000);
y_monte_cov400000 = cov(y_monte400000);  %reshape

% Store percentage mean error values in a matrix
Mean_per(:,i) = [ y_monte_mean300;    y_monte_mean25000;...
    y_monte_mean100000;   y_monte_mean400000; trig_mean_GenUT; trig_mean_UT];
Mean_per(:,i) = abs(Mean_per(:,i)- trig_mean_True)*100/trig_mean_True;
% Store percentage variance error values in a matrix
Cov_per(:,i) = [ y_monte_cov300;   y_monte_cov25000;...
    y_monte_cov100000;   y_monte_cov400000;  trig_cov_GenUT; trig_cov_UT];
Cov_per(:,i) = abs(Cov_per(:,i)- trig_cov_True)*100/trig_cov_True;
end
Mean_per = reshape(Mean_per,[6*num_iter,1]);
Cov_per = reshape(Cov_per,[6*num_iter,1]);
Method = repmat(Method,num_iter,1);
Method = char(Method); % Convert to character
% Generate box plot of the results
figure(1)
subplot(2,1,1)
boxplot(Mean_per,Method)
ylabel('% error of mean')

subplot(2,1,2)
boxplot(Cov_per,Method)
ylabel('% error of covariance')
xlabel('Method')
sgtitle('Weibull') 
 
%-------------------------------------------------------------------------%
%*********************  Poisson Random Variable *************************%
%-------------------------------------------------------------------------%

%\***** Poisson(lambda)         lambda = 0.1 ********\%
lambda = 0.1;      % rate parameter
[Poiss_mean, Poiss_second, Poiss_third, Poiss_fourth, Poiss_Ex ] =  Poisson_moments(lambda);
% Evalute and store results of the different unscented transforms
Sample_statistics(Poiss_mean, Poiss_second, Poiss_third, Poiss_fourth); 

% Evaluate and store the true mean and variance
E_sin_theta = exp(lambda*(cos(1) - 1))*sin(lambda*sin(1));   % E [  sin(theta) ]
E_cos_2_theta = exp(lambda*(cos(2) - 1))*cos(lambda*sin(2)); % E [  cos(2*theta) ]
[trig_mean_True,trig_cov_True] = analyticMean(E_sin_theta, E_cos_2_theta);

% Evaluate and store Monte Carlo draws
num_iter = 20;
Method = ["300 Monte"; "25000 Monte"; "100000 Monte"; "400000 Monte"; "GenUT" ; "UT"];
Mean_per = zeros(6,num_iter);   % Preallocate for storage
Cov_per = zeros(6,num_iter);    % Preallocate for storage
for i = 1:num_iter
% 300 Monte Carlo Draws
x = poissrnd(lambda*ones(1, 300));
y_monte300 = trigonometricTransform(x);
y_monte_mean300 = mean(y_monte300);
y_monte_cov300 = cov(y_monte300);

% 25000 Monte Carlo Draws
x = poissrnd(lambda*ones(1, 25000));
y_monte25000 = trigonometricTransform(x);
y_monte_mean25000 = mean(y_monte25000);
y_monte_cov25000 = cov(y_monte25000);

% 100000 Monte Carlo Draws
x = poissrnd(lambda*ones(1, 100000));
y_monte100000 = trigonometricTransform(x);
y_monte_mean100000 = mean(y_monte100000);
y_monte_cov100000 = cov(y_monte100000);

% 400000 Monte Carlo Draws
x = poissrnd(lambda*ones(1, 400000));
y_monte400000 = trigonometricTransform(x);
y_monte_mean400000 = mean(y_monte400000);
y_monte_cov400000 = cov(y_monte400000);  %reshape

% Store percentage mean error values in a matrix
Mean_per(:,i) = [ y_monte_mean300;    y_monte_mean25000;...
    y_monte_mean100000;   y_monte_mean400000; trig_mean_GenUT; trig_mean_UT];
Mean_per(:,i) = abs(Mean_per(:,i)- trig_mean_True)*100/trig_mean_True;
% Store percentage variance error values in a matrix
Cov_per(:,i) = [ y_monte_cov300;   y_monte_cov25000;...
    y_monte_cov100000;   y_monte_cov400000;  trig_cov_GenUT; trig_cov_UT];
Cov_per(:,i) = abs(Cov_per(:,i)- trig_cov_True)*100/trig_cov_True;
end
Mean_per = reshape(Mean_per,[6*num_iter,1]);
Cov_per = reshape(Cov_per,[6*num_iter,1]);
Method = repmat(Method,num_iter,1);
Method = char(Method); % Convert to character
% Generate box plot of the results
figure(2)
subplot(2,1,1)
boxplot(Mean_per,Method)
ylabel('% error of mean')

subplot(2,1,2)
boxplot(Cov_per,Method)
ylabel('% error of covariance')
xlabel('Method')
sgtitle('Poisson') 

 


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
    trig_cov_GenUT trig_cov_HOSPUT trig_cov_UT 
% Generate the ensembles and weights for the different transforms
[sigma_UT, weights_UT] = unscentedEnsemble(x_mu, x_cov, sqrt(3)); % scaled UT
[sigma_GenUT, weights_GenUT] = GenUT_Ensemble(x_mu, x_cov, x_skew, x_kurt); % GenUT
[sigma_HOSPUT, weights_HOSPUT] = HOSP(x_mu, x_cov, x_skew, x_kurt); % HOSPUT

% Evaluate and store scaled UT mean and variance
[trig_mean_UT,trig_cov_UT,~, ~] = Evaluate_sample_statistics(trigonometricTransform(sigma_UT),weights_UT);

% Evaluate and store GenUT mean and variance
[trig_mean_GenUT,trig_cov_GenUT,~, ~] = Evaluate_sample_statistics(trigonometricTransform(sigma_GenUT),weights_GenUT);

% Evaluate and store HOSPUT mean and variance
[trig_mean_HOSPUT ,trig_cov_HOSPUT ,~, ~] = Evaluate_sample_statistics(trigonometricTransform(sigma_HOSPUT),weights_HOSPUT);
end

