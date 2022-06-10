clear; clc;
% Paper Title: A Generalized Unscented Transformation for Probability Distributions
%
% This code simulates the results for Case Study 3 in our paper
addpath('Distribution Moments');
addpath('Unscented Transforms');

%%
% Generate correlated Poisson random numbers using sum of independent
% Poisson random numbers
mu3 = 0.2; % Third mean
mu1 = 9.8; mu2 = 1.8;  % First and Second Mean

mu_x1 = mu1 + mu3;  % Mean of I
mu_x2 = mu2 + mu3;  % Mean of R
P_x1_x2 = [mu_x1   mu3;  mu3  mu_x2]; % Covariance showing correlation

rho = mu3/ sqrt( (mu1 + mu3)*(mu2 + mu3)  ); % Pearson Correlation coefficient


%%
I = 10; R = 2; beta = 1.5; N = 100; gamma = 0.3;

% Analytic raw moments of Poissin(I) 
E_I_2 = I^2 + I;   % E[I^2]
E_I_3 = I^3 + 3*I^2 + I; % E[I^3]
E_I_4 = I^4 + 6*I^3 + 7*I^2 + I; % E[I^4]

% Analytic raw moments of Poissin(R) 
E_R_2 = R^2 + R;   % E[R^2]
E_R_3 = R^3 + 3*R^2 + R; % E[I^3]
E_R_4 = R^4 + 6*R^3 + 7*R^2 + R; % E[I^4]

% Mean of x1*x2. That is E(X1*X2)
E_x1_x2 = mu_x1*mu_x2  + mu3;   % E(X1*X2)

% Analytic mean of reduced SIR model 
y_mean_true = [I; R] + [beta*(N*I - E_I_2 - E_x1_x2)/N ;  gamma*I ];

%% Evaluate accuracy in propagating means and covariances
n = 2;                  % State dimension
mu = [I; R];            % Mean
P = P_x1_x2;           % Covariance
skewX = mu;             % Diagonal component of skewness
kurtX = 3*mu.^2 + mu;   % Diagonal component of kurtosis

% Monte Carlo Draws
nummonte = 10000000;
% Generate correlated Poisson random numbers using sum of independent
% Poisson random numbers
Y1 = poissrnd(mu1.*ones(1, nummonte)); Y2 = poissrnd(mu2.*ones(1, nummonte));  Y3 = poissrnd(mu3.*ones(1, nummonte));
X1 = Y1 + Y3;   X2 = Y2 + Y3; % New Poisson random numbers
x = [X1; X2];
y_monte100000 = nonlinTrans(x);
y_monte_mean100000 = mean(y_monte100000,2);
% covariance of reduced SIR model using 10^7 Monte Carlo draws
y_monte_cov100000 = cov(y_monte100000');
y_cov_true = y_monte_cov100000; % Use Monte Carlo as true covariance

% Standard UT Draws
[sigma_UT, weights_UT] = unscentedEnsemble(mu, P, sqrt(3));
y_UT = nonlinTrans(sigma_UT);
[y_UT_mean,y_UT_cov,~, ~] = Evaluate_sample_statistics(y_UT,weights_UT);

% GenUT approximation
[sigma_GUT, weights_GUT] = GenUT_Ensemble(mu, P, skewX, kurtX);
y_GUT = nonlinTrans(sigma_GUT);
[y_GUT_mean,y_GUT_cov,~, ~] = Evaluate_sample_statistics(y_GUT,weights_GUT);

% General HOSPUT approximation
[sigma_HOSPUT, weights_HOSPUT] = HOSP(mu, P, skewX, kurtX);
y_HOSPUT = nonlinTrans(sigma_HOSPUT);
[y_HOSPUT_mean,y_HOSPUT_cov,~, ~] = Evaluate_sample_statistics(y_HOSPUT,weights_HOSPUT);

% Calculate and display the mean errors
TotMean = [y_mean_true  y_GUT_mean  y_UT_mean  y_HOSPUT_mean  y_monte_mean100000];  
TotMean = abs(TotMean- y_mean_true)*100./y_mean_true; % Percentage error
disp('GenUT mean % error =');   disp(TotMean(:,2)); 
disp('UT mean % error =');      disp(TotMean(:,3));
disp('HOSPUT mean % error =');  disp(TotMean(:,4));
disp('Monte Carlo mean % error =');  disp(TotMean(:,5));

% Calculate and display covariance errors
TotCov = [y_cov_true    y_GUT_cov   y_UT_cov   y_HOSPUT_cov];
Acc_y_GUT_cov = abs(y_GUT_cov - y_cov_true)*100./y_cov_true; % Percentage error
Acc_y_UT_cov = abs(y_UT_cov - y_cov_true)*100./y_cov_true; % Percentage error
Acc_y_HOSPUT_cov = abs(y_HOSPUT_cov - y_cov_true)*100./y_cov_true; % Percentage error
Acc_y_monte_cov = abs(y_monte_cov100000 - y_cov_true)*100./y_cov_true; % Percentage error 
disp('GenUT covariance % error =');     disp(Acc_y_GUT_cov);
disp('UT covariance % error =');        disp(Acc_y_UT_cov);
disp('HOSPUT covariance % error =');    disp(Acc_y_HOSPUT_cov);
disp('Monte Carlo covariance % error =');    disp(Acc_y_monte_cov);

%%  Nonlinear function
function y = nonlinTrans(x)
I = 10; R = 2; beta = 1.5; N = 100; gamma = 0.3;
y = [I; R] + [beta*(N-x(1,:) - x(2,:)).*x(1,:)/N ;  gamma*x(1,:) ];
end
