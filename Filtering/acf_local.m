function r = acf_local(y, numLags)
%ACF_LOCAL Sample autocorrelation function.
%
%   r = acf_local(y)          lags 0..min(20, numel(y)-1)
%   r = acf_local(y, numLags) lags 0..numLags
%
% Drop-in replacement for the Econometrics Toolbox autocorr(), so the
% replication package runs on a base MATLAB install. Same biased estimator
% autocorr uses (the 1/N cancels between numerator and denominator), same
% default lag count, same orientation: r(1) == 1 and r(2) is the lag-1
% autocorrelation.

y = y(:);
N = numel(y);

if nargin < 2 || isempty(numLags)
    numLags = min(20, N - 1);
end

y = y - mean(y);
r = zeros(numLags + 1, 1);
for k = 0:numLags
    r(k + 1) = y(1:N - k)' * y(1 + k:N);
end
r = r / r(1);
end
