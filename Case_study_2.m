clear; clc;
% Paper Title: A Generalized Unscented Transformation for Probability Distributions
%
% This code simulates the results for Case Study 2 in our paper
addpath('Distribution Moments');
addpath('Unscented Transforms');

%%
% Poisson(0.1)
mu_Poiss = 0.1;               % Mean Poisson random var
P_poiss = mu_Poiss;           % Covariance of Poissom randon var
skewX_poiss = mu_Poiss;       % Diagonal component of skewness
kurtX_poiss = 3*mu_Poiss.^2 + mu_Poiss;   % Diagonal component of kurtosis

% Rayleigh (1)
sigma = 1;                      
mu_Ray = sigma*sqrt(pi/2);      % Mean of Rayleigh random var
P_Ray = (2 - pi/2)*sigma^2;     % Covariance of Rayleigh random var
skewX_Ray = sigma^3*sqrt(pi/2)*(pi - 3);       % Diagonal component of skewness
kurtX_Ray = -(sigma^4*(3*pi^2 - 32))/4;        % Diagonal component of kurtosis

% Combined mean, variance, skewness, and kurtosis
mu = [mu_Poiss; mu_Ray];            % Mean
P = diag([P_poiss  P_Ray]);         % Covariance
skewX = [skewX_poiss; skewX_Ray];   % Diagonal component of skewness
kurtX = [kurtX_poiss; kurtX_Ray];   % Diagonal component of kurtosis
n = length(mu);                     % Dimension of random vector

% Monte Carlo Draws
x1 = poissrnd(mu_Poiss.*ones(1, 10000000));
x2 = raylrnd(  sigma.*ones(1, 10000000)   );
x = [x1;x2];
y_monte100000 = nonlinTrans(x);
y_mean_true = mean(y_monte100000,2);
y_cov_true = cov(y_monte100000');

% Scaled UT Draws
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
TotMean = [y_mean_true  y_GUT_mean  y_UT_mean  y_HOSPUT_mean];  
TotMean = abs(TotMean- y_mean_true)*100./y_mean_true; % Percentage error
disp('GenUT mean % error =');   disp(TotMean(:,2)); 
disp('UT mean % error =');      disp(TotMean(:,3));
disp('HOSPUT mean % error =');  disp(TotMean(:,4));

% Calculate and display covariance errors
TotCov = [y_cov_true    y_GUT_cov   y_UT_cov   y_HOSPUT_cov];
Acc_y_GUT_cov = abs(y_GUT_cov - y_cov_true)*100./y_cov_true; % Percentage error
Acc_y_UT_cov = abs(y_UT_cov - y_cov_true)*100./y_cov_true; % Percentage error
Acc_y_HOSPUT_cov = abs(y_HOSPUT_cov - y_cov_true)*100./y_cov_true; % Percentage error
disp('GenUT covariance % error =');     disp(Acc_y_GUT_cov);
disp('UT covariance % error =');        disp(Acc_y_UT_cov);
disp('HOSPUT covariance % error =');    disp(Acc_y_HOSPUT_cov);


%%  Nonlinear function
function y = nonlinTrans(x)
y = [sin(x(1,:).*x(2,:)) ; cos(x(1,:).*x(2,:))];
end


 