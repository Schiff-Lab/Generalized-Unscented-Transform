% Paper Title: A Generalized Unscented Transformation for Probability Distributions
%
function [x,w] = unscentedEnsemble(mu,C,scale)
% Generates ensembles for the Scaled Unscented Transform
%
% Build an ensemble x and weights w that match the mean mu and covariance
% C with scaling parameter "scale" which is alpha*sqrt(n+kappa)
%
%
% Inputs
%
%   mu          -   Mean of random vector
%   C           -   Covariance
%   scale       -   Scaling parameter
%
% OUTPUTS
%
%   x           -   Matrix of sigma points
%   w           -   Vector of weights corresponding to each sigma point

    n=length(mu);
    if (nargin < 3)  scale = sqrt(n+1); end

    [U,S,~]=svd(C);
    
    U = U*diag(sqrt(diag(S)))*U';   %%% matrix square root of C
    %scale = alpha*sqrt(n+kappa);
    
    %%% Form the scaled unscented ensemble %%%
    x = zeros(n,2*n+1);
    x(:,1) = mu;
    x(:,(2:n+1)) = repmat(mu,1,n)+scale*U;
    x(:,(n+2:2*n+1)) = repmat(mu,1,n)-scale*U;
    
    %%% Determine scaled ensemble weights %%%
    w = [scale^2-n ones(1,2*n)/2]/scale^2;
    w = w.';
    
end

