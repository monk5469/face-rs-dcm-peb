function step07_parameter_recovery_mc(cfg)
%STEP07_PARAMETER_RECOVERY_MC Monte Carlo recovery of B-matrix parameters.

if nargin < 1 || isempty(cfg)
    cfg = default_config();
end

in_gcm = fullfile(cfg.output_dir, 'b_matrix_peb', 'GCM_peb_fitB.mat');
in_bma = fullfile(cfg.output_dir, 'b_matrix_peb', 'BMA_B.mat');
if ~exist(in_gcm, 'file') || ~exist(in_bma, 'file')
    error('B-matrix PEB outputs not found. Run step05_estimate_b_matrix_peb first.');
end

G = load(in_gcm, 'GCM_B');
B = load(in_bma, 'BMA_B');
GCM = G.GCM_B;
BMA = B.BMA_B;

sim_cfg = cfg;
sim_cfg.use_parfor = false;

snrs = cfg.simulation.snr_parameter_recovery;
n_rep = cfg.simulation.n_rep;
bad_subj = cfg.simulation.bad_subject_index;
figure_flag = cfg.simulation.figure_flag;

if cfg.use_parfor
    start_parallel_pool(cfg.parallel_workers);
end

out_dir = fullfile(cfg.output_dir, 'simulation_parameter_recovery');
ensure_dir(out_dir);

results_include = struct();
results_exclude = struct();

dq = parallel.pool.DataQueue;
total_iter = n_rep * numel(snrs);
clear update_progress_parameter_recovery
afterEach(dq, @(~) update_progress_parameter_recovery(total_iter));

for si = 1:numel(snrs)
    snr = snrs(si);
    r_include = nan(n_rep, 1);
    r_exclude = nan(n_rep, 1);
    mse_include = nan(n_rep, 1);
    mse_exclude = nan(n_rep, 1);

    parfor rep = 1:n_rep
        send(dq, rep);
        rng(cfg.simulation.rng_seed_parameter_recovery + rep + si * 1000);

        GCMsim = cell(size(GCM, 1), 1);
        for s = 1:size(GCM, 1)
            DCMsim = prepare_simulation_dcm(GCM{s}, BMA, []);

            if s == bad_subj
                DCMsim.Ep.A = DCMsim.Ep.A * cfg.simulation.bad_subject_scale;
                DCMsim.Ep.B = DCMsim.Ep.B * cfg.simulation.bad_subject_scale;
                DCMsim.Ep.C = DCMsim.Ep.C * cfg.simulation.bad_subject_scale;
            end

            [~, ~, DCMsim_generate] = spm_dcm_generate(DCMsim, snr, figure_flag);
            DCMsim_generate = apply_dcm_options(DCMsim_generate, cfg);
            GCMsim{s} = DCMsim_generate;
        end

        M = [];
        M.X = ones(size(GCMsim, 1), 1);
        M.Q = cfg.peb_precision;

        [~, PEBsim_inc] = fit_peb_model(GCMsim, M, cfg.peb_fields_b, sim_cfg);
        simBMA_inc = spm_dcm_peb_bmc(PEBsim_inc);

        orig_b = full(BMA.Ep(:));
        rec_b_inc = full(simBMA_inc.Ep(:));
        r_include(rep) = corr(orig_b, rec_b_inc);
        mse_include(rep) = mean((orig_b - rec_b_inc).^2);

        GCMsim_exc = GCMsim;
        GCMsim_exc(bad_subj) = [];
        M.X = ones(size(GCMsim_exc, 1), 1);

        [~, PEBsim_exc] = fit_peb_model(GCMsim_exc, M, cfg.peb_fields_b, sim_cfg);
        simBMA_exc = spm_dcm_peb_bmc(PEBsim_exc);
        rec_b_exc = full(simBMA_exc.Ep(:));
        r_exclude(rep) = corr(orig_b, rec_b_exc);
        mse_exclude(rep) = mean((orig_b - rec_b_exc).^2);
    end

    results_include(si).SNR = snr;
    results_include(si).r_mean = mean(r_include, 'omitnan');
    results_include(si).r_ci = prctile(r_include, [2.5 97.5]);
    results_include(si).r_all = r_include;
    results_include(si).mse_mean = mean(mse_include, 'omitnan');
    results_include(si).mse_all = mse_include;

    results_exclude(si).SNR = snr;
    results_exclude(si).r_mean = mean(r_exclude, 'omitnan');
    results_exclude(si).r_ci = prctile(r_exclude, [2.5 97.5]);
    results_exclude(si).r_all = r_exclude;
    results_exclude(si).mse_mean = mean(mse_exclude, 'omitnan');
    results_exclude(si).mse_all = mse_exclude;
end

save(fullfile(out_dir, 'parameter_recovery_mc_results.mat'), ...
    'results_include', 'results_exclude', 'cfg', '-v7.3');
end

function DCMsim = prepare_simulation_dcm(DCM, BMA, replacement_b)
keep_fields = {'Y', 'U', 'n', 'v', 'TE', 'delays', 'a', 'b', 'c', 'Ep'};
DCMsim = rmfield(DCM, setdiff(fieldnames(DCM), keep_fields));

if isempty(replacement_b)
    b_vec = full(BMA.Ep(:));
else
    b_vec = replacement_b(:);
end

prefix = zeros(DCMsim.n * DCMsim.n, 1);
b_full = [prefix; b_vec];
DCMsim.Ep.B = reshape(b_full, DCMsim.n, DCMsim.n, numel(b_full) / (DCMsim.n * DCMsim.n));
end

function update_progress_parameter_recovery(total)
persistent count;
if isempty(count)
    count = 0;
end
count = count + 1;
fprintf('Parameter recovery progress: %d / %d (%.2f%%)\n', count, total, 100 * count / total);
end
