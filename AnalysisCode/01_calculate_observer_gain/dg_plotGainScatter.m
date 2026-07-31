function dg_plotGainScatter(results)
% dg_plotGainScatter - scatter mean/median gain, moving vs stationary
%
%   dg_plotGainScatter(results)
%
% results is the table returned by DG_COMPAREGAINROI, with columns
% subject, meanGain_mov, medianGain_mov, meanGain_stat, medianGain_stat.
% Makes two scatter plots -- one for the mean, one for the median -- each
% with stationary on the x-axis, moving on the y-axis, a unity line, and
% each point labeled by subject.
%
% See also DG_COMPAREGAINROI

plotOne(results.meanGain_stat, results.meanGain_mov, results.subject, 'Mean gain');
plotOne(results.medianGain_stat, results.medianGain_mov, results.subject, 'Median gain');
return


% ------------------------------------------------------------------------
function plotOne(x, y, labels, ttl)
figure('Color', 'w');
scatter(x, y, 60, 'filled');
hold on;

lims = [0, max([x; y]) * 1.1];
plot(lims, lims, 'k--');
xlim(lims); ylim(lims);
axis square;

text(x + 0.015*diff(lims), y, labels, 'FontSize', 8, 'Interpreter', 'none');

xlabel('stationary gain (%BOLD)');
ylabel('moving gain (%BOLD)');
title(sprintf('%s: moving vs stationary (V1, vexpl>0.1, 4-8 deg)', ttl));
grid on;
return
