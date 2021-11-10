% Paper Title: A Generalized Unscented Transformation for Probability Distributions
%
function [x,weights] = HOSP(mu,P,x_skew, x_kurt)
% From the paper
%       A New Unscented Kalman Filter with Higher Order Moment-Matching
%
% Inputs
%
%   mu          -   Mean of random vector
%   P           -   Covariance
%   x_skew      -   Vector of diagonal components of the skewness tensor
%   x_kurt      -   Vector of diagonal components of the kurtosis tensor
%   lb          -   Lower bound of the state
%   ub          -   Upper bound of the state
%
% OUTPUTS
%
%   x           -   Matrix of sigma points
%   weights     -   Vector of weights corresponding to each sigma point


% Evaluate the matrix square root via singular value decomposition
[U1,S,~]=svd(P);
C = U1*diag(sqrt(diag(S)))*U1';   %%% matrix square root of C
N = length(mu);

phi1 = sum(x_skew)/ ( sqrt(N)*sum(sum(C.^3)) );
phi2 = sum(x_kurt)/ ( N*sum(sum(C.^4)) );

alpha = 0.5*( phi1 + sqrt(4*phi2 - 3*phi1^2)  ); 
beta = 0.5*( -phi1 + sqrt(4*phi2 - 3*phi1^2)  ) ;

% Calculate parameters
W = zeros(1,2*N);
n = N;
for i = 1:n
   W(i) = 1/(alpha*(alpha + beta)*N);
   W(n+i) = 1/(beta*(alpha + beta)*N);
end

% Generate sigma points
x0 = mu;
x1 = mu + alpha*sqrt(N)*C ;
x2 = mu - beta*sqrt(N)*C ;

x = [x0  x1  x2 ];

w0 = 1 - sum(W);
weights = [w0 W]';
end

