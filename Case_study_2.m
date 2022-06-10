clear; clc;
% Paper Title: A Generalized Unscented Transformation for Probability Distributions
%
% This code simulates the results for Case Study 2 in our paper
addpath('Distribution Moments');
addpath('Unscented Transforms');

%%
alpha1 = 1;     % shape parameter 1
alpha2 = 2;     % shape parameter 2
beeta = 0.3;      % scale parameter

alpha_x1 = alpha1;     % =  
alpha_x2 = alpha1 + alpha2;     % =   
E_x1_x2 =  alpha_x1*(alpha_x2 +1)*beeta^2; % E[ x1*x2 ]

% Gamma(1,0.3)
mu_x1 = alpha_x1*beeta;            % Mean Gamma random var
P_x1 = alpha_x1*beeta^2;           % Covariance of Gamma randon var
skew_x1 = 2*alpha_x1*beeta^3;      % Diagonal component of skewness
kurt_x1 =  3*alpha_x1*beeta^4*(alpha_x1 + 2);   % Diagonal component of kurtosis

% Gamma(3,0.3)                
mu_x2 = alpha_x2*beeta;            % Mean of Gamma random var
P_x2 = alpha_x2*beeta^2;           % Covariance of Gamma random var
skew_x2 = 2*alpha_x2*beeta^3;      % Diagonal component of skewness
kurt_x2 = 3*alpha_x2*beeta^4*(alpha_x2 + 2);        % Diagonal component of kurtosis

cross_cov =  E_x1_x2 - mu_x1*mu_x2; % Off-diagonal term of covariance matrix
% Combined mean, variance, skewness, and kurtosis
mu = [mu_x1; mu_x2];            % Mean
P =  [P_x1   cross_cov;  cross_cov  P_x2];         % Covariance Matrix
skewX = [skew_x1; skew_x2];     % Diagonal component of skewness
kurtX = [kurt_x1; kurt_x2];     % Diagonal component of kurtosis
n = length(mu);                 % Dimension of random vector

% Create correlated Gamma random vector via Monte Carlo Draws
nummonte = 10^7;
Temp1 = gamrnd(alpha1,beeta,1,nummonte); Temp2 = gamrnd(alpha2,beeta,1,nummonte);   
X1 = Temp1;   X2 = Temp1 + Temp2;   % New correlated Gamma random numbers
                                    % Correlation achieved via Temp1
                                            
x = [X1; X2];
y_monte100000 = nonlinTrans(x);
% y_mean_true = mean(y_monte100000,2);
E_sin_theta = sin(alpha_x1*atan(beeta)) / ( sqrt(1+ beeta^2) )^alpha_x1; % E [  sin(theta) ]
E_cos_theta = cos(-alpha_x2*atan(beeta)) / ( sqrt(1+ beeta^2) )^alpha_x2;  % E [  cos(theta) ]
y_mean_true = [E_sin_theta;E_cos_theta];
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
y = [sin(x(1,:)) ; cos(x(2,:))];
end
