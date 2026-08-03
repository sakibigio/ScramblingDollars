R = load('data/robustness_iota.mat'); r = R.rob_iota;
base = log(r(1).series.sigma_us_t);
fprintf('iota    lambda*   corr(log sig, base)   mean level shift (log)\n');
for i = 1:numel(r)
    ls = log(r(i).series.sigma_us_t);
    fprintf('%.4f  %.3f     %.5f              %+.3f\n', r(i).iota, r(i).lambda, corr(ls, base), mean(ls - base));
end
P = load('data/robustness_ploss.mat'); q = P.rob_ploss;
pb = log(q(6).series.sigma_us_t);
c = arrayfun(@(k) corr(log(q(k).series.sigma_us_t), pb), 1:numel(q));
fprintf('ploss sweep: min corr with baseline sigma = %.5f | lambda* range [%.3f, %.3f]\n', ...
    min(c), min([q.lambda]), max([q.lambda]));
