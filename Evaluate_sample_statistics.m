% Paper Title: A Generalized Unscented Transformation for Probability Distributions
%
function [sigma_mean,sigma_cov,sigma_skew, sigma_kurt] = Evaluate_sample_statistics(sigma_points,weights)
% This function is a tool that is used to evaluate the sample statistics 
% of sigma points of an unscented transform

% INPUTS
%
%   sigma_points        -   Matrix of sigma points
%   weights             -   Vector of weights corresponding to each sigma point
%
% OUTPUTS
%
%   sigma_mean          -   Sample mean of sigma points
%   sigma_cov           -   Sample covariance matrix of sigma points
%   sigma_skew          -   Sample diagonal component of skewness tensor
%   sigma_kurt          -   Sample diagonal component of kurtosis tensor

n = size(sigma_points,1); 
weights = weights.';        % Convert to row vector
% Mean
sigma_mean = sum(repmat(weights,n,1).*sigma_points,2); 
% Covariance
Z = (sigma_points - sigma_mean); 
sigma_cov = Z*diag(weights)*Z.';  
% Diagonal skewness
sigma_skew = sum(repmat(weights,n,1).*Z.^3,2);
% Diagonal kurtosis
sigma_kurt = sum(repmat(weights,n,1).*Z.^4,2);
end

