clear; clc;
% Paper Title: A Generalized Unscented Transformation for Probability Distributions
%
% This code simulates the results for Case Study 3 in our paper
addpath('Distribution Moments');
addpath('Unscented Transforms');

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

% Analytic mean of reduced SIR model 
y_mean_true = [I; R] + [beta*(N*I - E_I_2 - I*R)/N ;  gamma*I ];
% Analytic covariance of reduced SIR model 
cov_true_11 = I^2*beta^2 + E_I_2*beta^2 - (2*I^2*beta^2)/N + (I^2*beta^2)/N^2 - (2*I^3*beta^2)/N + (2*I^3*beta^2)/N^2 + (I^4*beta^2)/N^2 - (2*E_I_3*beta^2)/N + (E_I_4*beta^2)/N^2 - 2*I*I*beta^2 + (I^2*R^2*beta^2)/N^2 - (2*I^2*E_I_2*beta^2)/N^2 + (E_I_2*E_R_2*beta^2)/N^2 + (2*I*I*beta^2)/N - (2*I^2*R*beta^2)/N + (2*I^2*R*beta^2)/N^2 + (2*I^3*R*beta^2)/N^2 + (2*I*E_I_2*beta^2)/N + (2*I^2*I*beta^2)/N - (2*I*E_I_2*beta^2)/N^2 - (2*E_I_2*R*beta^2)/N + (2*E_I_3*R*beta^2)/N^2 - (2*I*R*E_I_2*beta^2)/N^2 - (2*I^2*I*R*beta^2)/N^2 + (2*I*R*I*beta^2)/N + (2*I*I*R*beta^2)/N - (2*I*I*R*beta^2)/N^2 - (2*I*R*I*R*beta^2)/N^2;
cov_true_22 = gamma^2*E_I_2 - 2*I*gamma^2*I + I^2*gamma^2;
cov_true_12 = I^2*beta*gamma + E_I_2*beta*gamma - (I^2*beta*gamma)/N - (I^3*beta*gamma)/N - (E_I_3*beta*gamma)/N - 2*I*I*beta*gamma + (I*I*beta*gamma)/N - (I^2*R*beta*gamma)/N + (I*E_I_2*beta*gamma)/N + (I^2*I*beta*gamma)/N - (E_I_2*R*beta*gamma)/N + (I*R*I*beta*gamma)/N + (I*I*R*beta*gamma)/N;
cov_true_21 = cov_true_12;
y_cov_true = [cov_true_11  cov_true_12;  cov_true_21  cov_true_22]; 

%% Evaluate accuracy in propagating means and covariances
n = 2;                  % State dimension
mu = [I; R];            % Mean
P = diag(mu);           % Covariance
skewX = mu;             % Diagonal component of skewness
kurtX = 3*mu.^2 + mu;   % Diagonal component of kurtosis

% Monte Carlo Draws
x = poissrnd(mu.*ones(2, 100000));
y_monte100000 = nonlinTrans(x);
y_monte_mean100000 = mean(y_monte100000,2);
y_monte_cov100000 = cov(y_monte100000');

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


 