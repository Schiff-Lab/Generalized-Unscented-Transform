clear;clc;
% Paper Title: A Generalized Unscented Transformation for Probability Distributions
%
% This code simulates the  results for Examples 1 and 2 in our paper
addpath('Distribution Moments');
addpath('Unscented Transforms');

xbar = [1.5;1];     % Mean
CovX = diag(xbar);  % Covariance of x

n = numel(xbar);
SkewX = xbar;               % Diagonal component of skewness of x
KurtX = 3*xbar.^2 + xbar;   % Diagonal component of kurtosis of x

% Call GenUT sigma point algorithm without constraints
[Xsigma, Warray_UPT,s_array] = GenUT_Ensemble(xbar, CovX, SkewX, KurtX);
% Calculate the sample mean and covariance
[sigma_mean,sigma_cov,diag_skew, diag_kurt] = Evaluate_sample_statistics(Xsigma,Warray_UPT);

%% Constrained Sigma Points of GenUT
% 
% Call constraied GenUT sigma point algorithm
lb = zeros(n,1);
[Xsigma_con, Warray_UPT_con,s_array_con] = GenUT_Ensemble(xbar, CovX, SkewX, KurtX, lb);
% Calculate the constrained mean and covariance
[sigma_mean_con,sigma_cov_con,con_skew, con_kurt] = Evaluate_sample_statistics(Xsigma_con,Warray_UPT_con);

%%
% Truncated version of sigma points
XsigmaT = max(Xsigma, lb);
% Calculate the truncated mean and covariance
[sigma_meanT,sigma_covT,T_skew, T_kurt] = Evaluate_sample_statistics(XsigmaT,Warray_UPT);

%%    Plot the sigma points
figure
scatter(Xsigma(1,:), Xsigma(2,:),100,'r', 'd','LineWidth',2);hold on;
scatter(XsigmaT(1,:), XsigmaT(2,:),100,'g', 'x','LineWidth',2);hold on;
scatter(Xsigma_con(1,:), Xsigma_con(2,:),100,'b', 'c','LineWidth',2)
xlabel('x_1'); ylabel('x_2'); grid on;
legend('Unconstrained', 'Truncated', 'Constrained')
xlim([-0.5 4.5])               % Set x-axis limits
ylim([-0.5 3.5])               % Set x-axis limits


% Plot the covariances against the true one
a = sqrt(CovX(1,1))/2; % horizontal radius
b = sqrt(CovX(2,2))/2; % vertical radius
x0 = xbar(1); y0 = xbar(2);   t = -pi:0.01:pi;
x1_true = x0+a*cos(t); x2_true = y0+b*sin(t);  % True covariance

a = sqrt(sigma_cov(1,1))/2; % horizontal radius
b = sqrt(sigma_cov(2,2))/2; % vertical radius
x0 = sigma_mean(1); y0 = sigma_mean(2);   t = -pi:0.01:pi;
x1_unc = x0+a*cos(t); x2_unc = y0+b*sin(t);  % Unconstrained covariance

a = sqrt(sigma_covT(1,1))/2; % horizontal radius
b = sqrt(sigma_covT(2,2))/2; % vertical radius
x0 = sigma_meanT(1); y0 = sigma_meanT(2);   t = -pi:0.01:pi;
x1_Trunc = x0+a*cos(t); x2_Trunc = y0+b*sin(t);  % Truncated covariance

a = sqrt(sigma_cov_con(1,1))/2; % horizontal radius
b = sqrt(sigma_cov_con(2,2))/2; % vertical radius
x0 = sigma_mean_con(1); y0 = sigma_mean_con(2);   t = -pi:0.01:pi;
x1_con = x0+a*cos(t); x2_con = y0+b*sin(t);  % Constrained covariance

figure
scatter(xbar(1), xbar(2), 100,'k', '+','LineWidth',3); hold on;
scatter(sigma_mean(1), sigma_mean(2), 100,'r', 'd','LineWidth',3); hold on;
scatter(sigma_meanT(1), sigma_meanT(2), 100,'g', 'x','LineWidth',3); hold on;
scatter(sigma_mean_con(1), sigma_mean_con(2), 100,'b', 'c','LineWidth',3); hold on;
plot(x1_true,x2_true, 'k', 'LineWidth',3); hold on;
plot(x1_unc,x2_unc, '-.r', 'LineWidth',3); hold on;
plot(x1_Trunc,x2_Trunc, ':g', 'LineWidth',3); hold on;
plot(x1_con,x2_con, '--b', 'LineWidth',3);
xlabel('x_1'); ylabel('x_2');
xlim([0.8 2.2])               % Set x-axis limits
ylim([0.4 1.6])               % Set x-axis limits
legend('True Mean', 'Unconstrained Mean','Truncated Mean', 'Constrained Mean',...
   'True Covariance', 'Unconstrained Covariance','Truncated Covariance', 'Constrained Covariance');
grid