function step04_loo_input_model_validation(cfg)
%STEP04_LOO_INPUT_MODEL_VALIDATION Leave-one-out validation for input BMC.

if nargin < 1 || isempty(cfg)
    cfg = default_config();
end

in_file = fullfile(cfg.output_dir, 'input_model_comparison', 'GCM_fit.mat');
if ~exist(in_file, 'file')
    error('Fitted GCM file not found: %s', in_file);
end

S = load(in_file, 'GCM', 'input_models');
GCM = S.GCM;
input_models = S.input_models;
F = extract_free_energy_matrix(GCM);
[n_sub, n_mod] = size(F);

[~, exp_r_full, xp_full, pxp_full, bor_full] = spm_BMS(F);
[~, winner_full] = max(pxp_full);

loo_post_pred = nan(n_sub, n_mod);
loo_exp_r = nan(n_sub, n_mod);
loo_winner_pred = nan(n_sub, 1);
loo_winner_true = nan(n_sub, 1);

for i = 1:n_sub
    F_train = F([1:i-1, i+1:end], :);
    [alpha_loo, exp_r_loo] = spm_BMS(F_train);
    loo_exp_r(i, :) = exp_r_loo;

    F_i = F(i, :);
    F_i_rel = F_i - max(F_i);
    likelihood = exp(F_i_rel);
    post_pred = alpha_loo .* likelihood;
    post_pred = post_pred / sum(post_pred);

    loo_post_pred(i, :) = post_pred;
    [~, loo_winner_pred(i)] = max(post_pred);
    [~, loo_winner_true(i)] = max(F_i);
end

correct = loo_winner_pred == loo_winner_true;
loo_accuracy_pct = mean(correct) * 100;

out_dir = fullfile(cfg.output_dir, 'input_model_comparison', 'loo_validation');
ensure_dir(out_dir);

subject_id = arrayfun(@(s) sprintf('sub-%02d', s), (1:n_sub)', 'UniformOutput', false);
true_model = model_names_from_index(input_models, loo_winner_true);
pred_model = model_names_from_index(input_models, loo_winner_pred);

T_loo = table(subject_id, loo_winner_true, loo_winner_pred, correct, ...
    loo_post_pred(:, 1), loo_post_pred(:, 2), loo_post_pred(:, 3), ...
    true_model, pred_model, ...
    'VariableNames', {'SubjectID', 'TrueModel_idx', 'PredModel_idx', 'Correct', ...
    'PostPred_M1', 'PostPred_M2', 'PostPred_M3', 'TrueModel', 'PredModel'});
writetable(T_loo, fullfile(out_dir, 'LOO_predictions.csv'));

T_summary = table(loo_accuracy_pct, sum(correct), n_sub, bor_full, ...
    pxp_full(winner_full), xp_full(winner_full), exp_r_full(winner_full), ...
    string(input_models{winner_full}.name), ...
    'VariableNames', {'LOO_Accuracy_pct', 'N_Correct', 'N_Total', ...
    'FullSample_BOR', 'FullSample_max_pxp', 'FullSample_max_xp', ...
    'FullSample_max_exp_r', 'WinningModel'});
writetable(T_summary, fullfile(out_dir, 'LOO_summary.csv'));

save(fullfile(out_dir, 'LOO_validation.mat'), 'F', 'loo_post_pred', ...
    'loo_exp_r', 'loo_winner_pred', 'loo_winner_true', 'correct', ...
    'loo_accuracy_pct', 'bor_full', 'pxp_full', 'xp_full', 'exp_r_full', ...
    'input_models', 'cfg');
end

function names = model_names_from_index(input_models, model_index)
names = strings(numel(model_index), 1);
for i = 1:numel(model_index)
    if isnan(model_index(i))
        names(i) = missing;
    else
        names(i) = string(input_models{model_index(i)}.name);
    end
end
end
