function step03_define_fit_input_models(cfg)
%STEP03_DEFINE_FIT_INPUT_MODELS Define and fit DCM input-model space.

if nargin < 1 || isempty(cfg)
    cfg = default_config();
end

paths = load_workspace_paths(cfg);
input_models = build_input_models();
voi_names = make_voi_names(cfg);
n_voi = numel(cfg.voi_indices);

out_root = fullfile(cfg.output_dir, 'input_model_comparison');
ensure_dir(out_root);

GCM = cell(paths.n_subjects * cfg.session_count, numel(input_models));
row_index = 0;

for s = 1:paths.n_subjects
    spm_file = fullfile(paths.target_dir{s}, 'bold', cfg.first_level_dir_name, 'SPM.mat');
    if ~exist(spm_file, 'file')
        error('SPM.mat not found for subject %d: %s', s, spm_file);
    end
    S = load(spm_file, 'SPM');
    SPM = S.SPM;

    for ses = 1:cfg.session_count
        row_index = row_index + 1;
        xY = repmat(struct(), 1, n_voi);

        for v = 1:n_voi
            roi_index = cfg.voi_indices(v);
            voi_file = fullfile(paths.target_dir{s}, 'bold', cfg.roi_dir_name, ...
                [voi_names{roi_index, ses} '.mat']);
            if ~exist(voi_file, 'file')
                error('VOI file not found for subject %d: %s', s, voi_file);
            end
            V = load(voi_file, 'xY');
            xY(v) = V.xY;
        end

        for m = 1:numel(input_models)
            DCM = build_dcm_from_spm(SPM, xY, input_models{m}, cfg, ses);

            model_dir = fullfile(paths.target_dir{s}, 'bold', cfg.dcm_dir_name, input_models{m}.name);
            ensure_dir(model_dir);
            dcm_file = fullfile(model_dir, sprintf('DCM_mod1_ses%d.mat', ses));
            save(dcm_file, 'DCM');
            GCM{row_index, m} = dcm_file;
        end
    end
end

save(fullfile(out_root, 'GCM_defined.mat'), 'GCM', 'input_models', 'cfg', '-v7.3');

if cfg.use_parfor
    start_parallel_pool(cfg.parallel_workers);
end

GCM = spm_dcm_load(GCM);
GCM = spm_dcm_fit(GCM, cfg.use_parfor);
save(fullfile(out_root, 'GCM_fit.mat'), 'GCM', 'input_models', 'cfg', '-v7.3');

F = extract_free_energy_matrix(GCM);
[alpha, exp_r, xp, pxp, bor] = spm_BMS(F);
save(fullfile(out_root, 'BMC_input_models.mat'), ...
    'F', 'alpha', 'exp_r', 'xp', 'pxp', 'bor', 'input_models', 'cfg');

M = [];
M.X = ones(size(GCM, 1), 1);
M.Q = cfg.peb_precision;
[PEB_C, RCM_C] = spm_dcm_peb(GCM(:, cfg.selected_input_model), M, cfg.peb_fields_c);
BMA_C = spm_dcm_peb_bmc(PEB_C);
save(fullfile(out_root, 'PEB_C_selected_input_model.mat'), ...
    'PEB_C', 'RCM_C', 'BMA_C', 'cfg', '-v7.3');
end
