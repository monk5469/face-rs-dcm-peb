function run_pipeline(cfg)
%RUN_PIPELINE Minimal reproducible entry point for the DCM input analysis.

if nargin < 1 || isempty(cfg)
    cfg = default_config();
end

if isfield(cfg, 'run_all') && cfg.run_all
    if isfield(cfg, 'run_preprocessing') && cfg.run_preprocessing
        step01_preprocess_fmri(cfg);
    end
    if isfield(cfg, 'run_glm') && cfg.run_glm
        step02_first_level_glm(cfg);
    end
    step03_define_fit_input_models(cfg);
    step04_loo_input_model_validation(cfg);
    step05_estimate_b_matrix_peb(cfg);
    step06_compute_variance_explained(cfg);
    step07_parameter_recovery_mc(cfg);
    step08_detection_sensitivity_power(cfg);
else
    step03_define_fit_input_models(cfg);
end
end
