function [X] = logistic_reflect(n, W, R, C, X0)
%% generate time series from coupled logistic map with reflecting borders
% n ... length of desired time series
% W ... dxd-dimensional matrix of coupling weights
% R ... dx1 dimensional parameter for logistic maps 
% C ... dxd dimensional correlation matrix for interal noise
% X0... 1xd-dimensional starting values (optional), default are uniformly dsitributed random values from  interval [0,1]
% 
% X ... nxd output-timeseries

%% evaluate input
% determine dimension d
if iscolumn(R) == 1
    R = R';
end
d = size(R, 2);

% check dimension of weights
if size(W, 1) ~= d || size(W, 2) ~= d 
    error('Error: weight matrix does not match dimension')
end

% set starting values if needed
if nargin < 5
    X0 = rand(1, d);
else
    if iscolumn(X0) == 1
        X0 = X0';
    end
    if size(X0, 2) ~= d || size(X0, 1) ~= 1
        error('Error: Starting values do not match dimension!')
    end
end

%% generate time-series
X = zeros(n, d);
noise = randn(n, d);
X(1, :) = X0;
for t = 2:n
    X_h = X(t-1, :).*(R.*(1 - X(t-1, :)) + X(t-1, :)*W') + noise(t-1, :)*C';
    X(t, :) = reflect(X_h);
end
end

function [X] = reflect(Xi)
% reflect values leaving the interval [0,1] back
if (sum(Xi >= 1) > 0) || (sum(Xi <= 0) > 0)
    % reflect negativ values
    Xi = -Xi.*(Xi<=0) + Xi.*(Xi>0);
    % reflect values > 1
    X = (1 - Xi + (Xi>=1)).*(Xi>=1) + Xi.*(Xi<1); 
else
    X = Xi;
end
end
