function step05_estimate_b_matrix_peb(cfg)
%STEP05_ESTIMATE_B_MATRIX_PEB Estimate B-matrix modulation with PEB/BMA.

if nargin < 1 || isempty(cfg)
    cfg = default_config();
end

in_file = fullfile(cfg.output_dir, 'input_model_comparison', 'GCM_fit.mat');
if ~exist(in_file, 'file')
    error('Fitted GCM file not found: %s', in_file);
end

S = load(in_file, 'GCM', 'input_models');
GCM = S.GCM(:, cfg.selected_input_model);
input_model = S.input_models{cfg.selected_input_model};
if ischar(GCM{1}) || isstring(GCM{1})
    GCM = spm_dcm_load(GCM);
end

out_dir = fullfile(cfg.output_dir, 'b_matrix_peb');
ensure_dir(out_dir);

M = [];
M.X = ones(size(GCM, 1), 1);
M.Q = cfg.peb_precision;

[GCM_B, PEB_B] = fit_peb_model(GCM, M, cfg.peb_fields_b, cfg);
BMA_B = spm_dcm_peb_bmc(PEB_B);
[FamP, family_table, family_results] = compute_b_family_probabilities( ...
    PEB_B, GCM_B, input_model);

save(fullfile(out_dir, 'GCM_peb_fitB.mat'), 'GCM_B', 'PEB_B', 'input_model', 'cfg', '-v7.3');
save(fullfile(out_dir, 'BMA_B.mat'), 'BMA_B', 'input_model', 'cfg', '-v7.3');
save(fullfile(out_dir, 'B_family_probabilities.mat'), ...
    'FamP', 'family_table', 'family_results', 'input_model', 'cfg', '-v7.3');
writetable(family_table, fullfile(out_dir, 'B_family_posterior_probabilities.csv'));

if isfield(BMA_B, 'Ep') && isfield(BMA_B, 'Pp')
    n_params = numel(full(BMA_B.Ep));
    parameter_index = (1:n_params)';
    posterior_expectation = full(BMA_B.Ep(:));
    posterior_probability = full(BMA_B.Pp(:));
    selected = posterior_probability >= cfg.bma_probability_threshold;

    T = table(parameter_index, posterior_expectation, posterior_probability, selected, ...
        'VariableNames', {'ParameterIndex', 'Ep_Hz', 'Pp', 'Selected_Pp_Threshold'});
    writetable(T, fullfile(out_dir, 'BMA_B_parameter_summary.csv'));
end
end

function [FamP, family_table, family_results] = compute_b_family_probabilities(PEB, GCM, input_model)
vals = [0 0; 0 1; 1 0; 1 1];
n_val = size(vals, 1);
n_models = n_val ^ 2;

family_names = { ...
    'Any self-connection modulation', ...
    'Any between-region modulation', ...
    'rOFA self-connection modulation', ...
    'rFFA self-connection modulation', ...
    'Feedforward rOFA_to_rFFA modulation', ...
    'Feedback rFFA_to_rOFA modulation'};

n_family = numel(family_names);
n_inputs = size(input_model.contrast, 1);
n_driver = numel(input_model.driver_rows);
n_modulatory = n_inputs - n_driver;

FamP = nan(n_family, n_modulatory);
family_results = cell(n_family, n_modulatory);

base_model = GCM{1, 1};
if isfield(base_model, 'M')
    base_model = rmfield(base_model, 'M');
end

row_condition = [];
row_input = [];
row_condition_name = {};
row_family = [];
row_family_name = {};
row_probability = [];

final_PF = [];
for con = n_driver + 1:n_inputs
    GCMs = cell(1, n_models);
    m = 1;
    GCMs{1, m} = base_model;

    family = nan(n_family, n_models);
    family(:, 1) = ones(n_family, 1) * 2;

    for self = 1:n_val
        for ofa_ffa = 1:n_val
            if m < n_models
                m = m + 1;
                GCMs{1, m} = base_model;

                GCMs{1, m}.b(1, 1, con) = vals(self, 1);
                GCMs{1, m}.b(2, 2, con) = vals(self, 2);
                GCMs{1, m}.b(2, 1, con) = vals(ofa_ffa, 1);
                GCMs{1, m}.b(1, 2, con) = vals(ofa_ffa, 2);

                family(1, m) = max(vals(self, :)) + 1;
                family(2, m) = max(vals(ofa_ffa, :)) + 1;
                family(3, m) = vals(self, 1) + 1;
                family(4, m) = vals(self, 2) + 1;
                family(5, m) = vals(ofa_ffa, 1) + 1;
                family(6, m) = vals(ofa_ffa, 2) + 1;
            end
        end
    end

    [BMA_family, BMR_family] = spm_dcm_peb_bmc(PEB, GCMs);

    for f = 1:n_family
        [~, fam_task, final_PF] = spm_dcm_peb_bmc_fam( ...
            BMA_family, BMR_family, family(f, :), 'ALL', final_PF);
        family_results{f, con - n_driver} = fam_task;
        FamP(f, con - n_driver) = fam_task.family.post(2);

        row_condition(end + 1, 1) = con - n_driver; %#ok<AGROW>
        row_input(end + 1, 1) = con; %#ok<AGROW>
        row_condition_name{end + 1, 1} = input_name(input_model, con); %#ok<AGROW>
        row_family(end + 1, 1) = f; %#ok<AGROW>
        row_family_name{end + 1, 1} = family_names{f}; %#ok<AGROW>
        row_probability(end + 1, 1) = FamP(f, con - n_driver); %#ok<AGROW>
    end
end

family_table = table(row_condition, row_input, row_condition_name, ...
    row_family, row_family_name, row_probability, ...
    'VariableNames', {'ModulatoryConditionIndex', 'DCMInputIndex', ...
    'ConditionName', 'FamilyIndex', 'FamilyName', 'PosteriorProbability'});
end

function name = input_name(input_model, con)
condition_labels = {'FF_Initial', 'FF_Immediate', 'FF_Delayed', ...
    'UF_Initial', 'UF_Immediate', 'UF_Delayed', ...
    'SF_Initial', 'SF_Immediate', 'SF_Delayed'};

parts = condition_labels(input_model.contrast(con, :) ~= 0);
name = strjoin(parts, '+');
end
