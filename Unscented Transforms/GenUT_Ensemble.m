% Paper Title: A Generalized Unscented Transformation for Probability Distributions
%
% This code generates sigma points using our Generalized Unscented Transform
function [x,weights,s] = GenUT_Ensemble(mu,P,x_skew, x_kurt, lb, ub)
% INPUTS
%
%   mu          -   Mean of random vector
%   P           -   Covariance matrix
%   x_skew      -   Vector of diagonal components of the skewness tensor
%   x_kurt      -   Vector of diagonal components of the kurtosis tensor
%   lb          -   Vector of lower bound of the state
%   ub          -   Vector of upper bound of the state
%
% OUTPUTS
%
%   x           -   Matrix of sigma points
%   weights     -   Vector of weights corresponding to each sigma point
%   s           -   Vector of proportionality constants used to generate
%                   the sigma points

%   For the bounds, we note that  lb < x  < ub

% Get the number of states
n = size(mu,1);

% Evaluate the matrix square root via singular value decomposition
[U1,S,~]=svd(P);
C = U1*diag(sqrt(diag(S)))*U1';   % C is the matrix square root of P

% Handle the arguments for skewness and kurtosis
if nargin < 3 || isempty(x_skew)    % If no diagonal component of skewness
    % is specified
    warning('No skewness specified: Gaussian skewness and kurtosis is assumed')
    x_skew = 0*mu;          % Assume gaussian skewness if not provided
    x_kurt = 3*ones(n,1);   % Assume gaussian diagonal kurtosis in turn
end
if nargin < 4 || isempty(x_kurt)    % If no diagonal component of kurtosis
    % is specified
    warning('No kurtosis specified: kurtosis is selected to satisfy skewness kurtosis relationship')
    % Manually ensure kurtosis skew relationship is satisfied
    x_kurt = (C.^4)*((C.^3)\x_skew).^2;       x_kurt = 1.1*x_kurt;
end

% Handle when specified kurtosis violates skewness kurtosis relationship
minkurt = (C.^4)*((C.^3)\x_skew).^2;
if (sum(x_kurt<minkurt))
    warning('Bad Human Error: Kurtosis does not correspond to a distribution')
    x_kurt(x_kurt<minkurt) = 1.001*minkurt(x_kurt<minkurt);
end

% Handle the arguments for lower bounds and upper bounds
if nargin < 5 || isempty(lb) % If lower bound is not specified
    % Manually set lower bound as -inf
    lb = -inf*ones(n,1);
end
if nargin < 6 || isempty(lb) % If upper bound is not specified
    % Manually set upper bound as inf
    ub = inf*ones(n,1);
end

% Calculate parameters u and v
u = 0.5*( -((C.^3)\x_skew) + ...
    sqrt( 4*((C.^4)\x_kurt) - 3*((C.^3)\x_skew).^2 )    );
v = u + (C.^3)\x_skew;

% Generate the sigma points
x0 = mu;
x1 = mu - C*diag(u);
x2 = mu + C*diag(v);

%% --------------- This section handles the constraints  --------------- %%
Flag_constrain = 0;     % Default flag to enforce constraint
% Check if mean violates constraints
if min(mu - lb)<0 || min(ub - mu)<0
    Flag_constrain = 1; % Set flag to avoid enforcing state constraints
    warning('Unable fo enforce constraints: one or more of the mean does not satisfy lb < mean < ub')
end

if Flag_constrain == 0
    theta = 0.9;    % Default value of user defined slack parameter
    
    %\     Ensure lower bound 'lb' is not violated
    Temp1 = [x1   x2] - lb;
    L1 = find(min(Temp1)<0);    % Find the location of sigma points that
                                % violate the lower bound
    Flag_calc = 0;      % Flag that determines if skewness can be matched
    for i = 1:length(L1)
        if L1(i) <= n
            % Recalculate 'u' to satisfy lower bound 'lb'
            u(L1(i)) =  theta*min(   abs(  (mu - lb)./ C(:,L1(i))   )   );
        else
            % Recalculate 'v' to satisfy lower bound 'lb'
            v(L1(i)-n) =   theta*min(   abs(  (lb - mu)./ C(:,L1(i)-n)  ));
            Flag_calc = 1;   % Set flag
        end
    end
    
    % Regenerate the sigma points
    x1 = mu - C*diag(u);
    x2 = mu + C*diag(v);
    
    %\     Ensure upper bound 'ub' is not violated
    Temp2 = ub  - [x1  x2];
    L2 = find(min(Temp2)<0);    % Find the location of sigma points that
                                % violate the upper bound
    for i = 1:length(L2)
        if L2(i) <= n
            % Recalculate 'u' to satisfy upper bound 'ub'
            u(L2(i)) =  theta*min(   abs(  (mu - ub)./  C(:,L2(i))   )   );
        else
            % Recalculate 'v' to satisfy upper bound 'ub'
            v(L2(i)-n) =   theta*min(   abs(  (ub - mu)./  C(:,L2(i)-n)   )   );
            Flag_calc = 1;   % Set flag
        end
        
    end
    
    if Flag_calc == 0
        % Now recalculate parameter 'v' to match diagonal componen of 
        % skewness tensor because it was 'v' was not previously redefined
        v = u + (C.^3)\x_skew; % only done of v was not redefined
    end
    
    % Regenerate the sigma points to reflect any change in 'u' or 'v'
    x1 = mu - C*diag(u);
    x2 = mu + C*diag(v);
end 

% Recalculate weights to reflect any change in 'u' or 'v'
w2 = ones(n,1)./ v ./ (u + v);
w1 = w2 .* v ./ u;

% Output sigma point values
x = [x0  x1  x2];
w0 = 1 - sum([w1; w2]);
weights = [w0 w1'  w2']';
s = [u; v];
end

