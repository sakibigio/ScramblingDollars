%% build_corridor_tbill.m — one-shot: build the DW - 1M T-bill corridor
% (base_mode='tbill' variant of build_iota_corridor_monthly, matching the
% Appendix H.1 specification) and save it to data/iota_corridor_tbill.mat,
% preserving the existing spliced-base file untouched.

if ~exist('data/iota_corridor_monthly.mat', 'file')
    error('expected existing data/iota_corridor_monthly.mat (spliced) to preserve');
end
copyfile('data/iota_corridor_monthly.mat', 'data/iota_corridor_monthly_SPLICED_BAK.mat');

base_mode = 'tbill';
run('build_iota_corridor_monthly.m');

movefile('data/iota_corridor_monthly.mat', 'data/iota_corridor_tbill.mat');
movefile('data/iota_corridor_monthly_SPLICED_BAK.mat', 'data/iota_corridor_monthly.mat');
fprintf('\n[build_corridor_tbill] wrote data/iota_corridor_tbill.mat; spliced file restored.\n');
