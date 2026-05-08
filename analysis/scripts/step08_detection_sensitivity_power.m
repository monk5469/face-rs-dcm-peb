function step08_detection_sensitivity_power(cfg)
%STEP08_DETECTION_SENSITIVITY_POWER Estimate posterior detection sensitivity.

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

snrs = cfg.simulation.snr_detection;
n_rep = cfg.simulation.n_rep;
bad_subj = cfg.simulation.bad_subject_index;
target_conns = cfg.simulation.target_connection_indices;
sigma_subj = cfg.simulation.subject_jitter_sd;
figure_flag = cfg.simulation.figure_flag;

if cfg.use_parfor
    start_parallel_pool(cfg.parallel_workers);
end

out_dir = fullfile(cfg.output_dir, 'simulation_detection_sensitivity');
ensure_dir(out_dir);

results = struct();
dq = parallel.pool.DataQueue;
total_iter = n_rep * numel(snrs);
clear update_progress_detection_sensitivity
afterEach(dq, @(~) update_progress_detection_sensitivity(total_iter));

for si = 1:numel(snrs)
    snr = snrs(si);

    for ci = 1:numel(target_conns)
        conn_idx = target_conns(ci);
        fitted_b = full(BMA.Ep(:));
        b_obs = fitted_b(conn_idx);
        multipliers = choose_effect_multipliers(b_obs);
        effect_grid = sign(b_obs) .* (multipliers * abs(b_obs));

        power_include = nan(numel(effect_grid), 1);
        power_exclude = nan(numel(effect_grid), 1);

        for ei = 1:numel(effect_grid)
            effect = effect_grid(ei);
            detected_include = false(n_rep, 1);
            detected_exclude = false(n_rep, 1);

            parfor rep = 1:n_rep
                send(dq, rep);
                rng(cfg.simulation.rng_seed_detection + si * 1e6 + ci * 1e4 + ei * 100 + rep);

                subj_jitter = normrnd(0, sigma_subj, [size(GCM, 1), 1]);
                GCMsim = cell(size(GCM, 1), 1);

                for s = 1:size(GCM, 1)
                    b_replaced = fitted_b;
                    b_replaced(conn_idx) = effect * (1 + subj_jitter(s));
                    DCMsim = prepare_detection_dcm(GCM{s}, b_replaced);

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
                detected_include(rep) = full(simBMA_inc.Pp(conn_idx)) > cfg.bma_probability_threshold;

                GCMsim_exc = GCMsim;
                GCMsim_exc(bad_subj) = [];
                M.X = ones(size(GCMsim_exc, 1), 1);
                [~, PEBsim_exc] = fit_peb_model(GCMsim_exc, M, cfg.peb_fields_b, sim_cfg);
                simBMA_exc = spm_dcm_peb_bmc(PEBsim_exc);
                detected_exclude(rep) = full(simBMA_exc.Pp(conn_idx)) > cfg.bma_probability_threshold;
            end

            power_include(ei) = mean(detected_include);
            power_exclude(ei) = mean(detected_exclude);
        end

        results(si, ci).SNR = snr;
        results(si, ci).conn = conn_idx;
        results(si, ci).conn_eff = b_obs;
        results(si, ci).multipliers = multipliers;
        results(si, ci).eff_grid = effect_grid;
        results(si, ci).power_include = power_include;
        results(si, ci).power_exclude = power_exclude;
    end
end

save(fullfile(out_dir, 'detection_sensitivity_results.mat'), 'results', 'cfg', '-v7.3');
end

function multipliers = choose_effect_multipliers(b_obs)
if abs(b_obs) < 0.3
    multipliers = [0.5 0.6 0.7 1.25 1.5 1.75 2 2.25 2.5 2.75 3];
elseif abs(b_obs) < 1
    multipliers = [0.5 0.52 0.55 0.58 0.6 0.7 1 1.25];
else
    multipliers = [0.2 0.23 0.25 0.28 0.3 0.5 0.7 1];
end
end

function DCMsim = prepare_detection_dcm(DCM, b_vec)
keep_fields = {'Y', 'U', 'n', 'v', 'TE', 'delays', 'a', 'b', 'c', 'Ep'};
DCMsim = rmfield(DCM, setdiff(fieldnames(DCM), keep_fields));
prefix = zeros(DCMsim.n * DCMsim.n, 1);
b_full = [prefix; b_vec(:)];
DCMsim.Ep.B = reshape(b_full, DCMsim.n, DCMsim.n, numel(b_full) / (DCMsim.n * DCMsim.n));
end

function update_progress_detection_sensitivity(total)
persistent count;
if isempty(count)
    count = 0;
end
count = count + 1;
fprintf('Detection sensitivity progress: %d / %d (%.2f%%)\n', count, total, 100 * count / total);
end
