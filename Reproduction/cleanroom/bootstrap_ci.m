function [ci, boots] = bootstrap_ci(x, nBoot, level, seed)
% BOOTSTRAP_CI  Bootstrap CI of the mean of a per-subject statistic.
%
%   [ci, boots] = bootstrap_ci(x, nBoot, level, seed)
%
% x     : 1 x nSubj vector (one statistic per subject).
% nBoot : number of resamples (default 1000).
% level : CI percent (default 95).
% seed  : rng seed for reproducibility (default 0).
%
% Resamples subjects with replacement; returns ci = [lo hi] percentile interval
% and the bootstrap distribution of the mean.

    if nargin < 2 || isempty(nBoot), nBoot = 1000; end
    if nargin < 3 || isempty(level), level = 95;   end
    if nargin < 4 || isempty(seed),  seed  = 0;    end

    rng(seed);
    x = x(:).';
    n = numel(x);
    boots = zeros(nBoot, 1);
    for b = 1:nBoot
        boots(b) = mean(x(randi(n, [1 n])));
    end
    a = (100 - level) / 2;
    ci = prctile(boots, [a, 100 - a]);
end
