# Standardized Analysis Scripts

This folder contains MATLAB scripts for the DCM-PEB analyses.

The scripts are written as functions with a shared `cfg` structure. Start by editing `config/default_config.m`, then run individual steps from MATLAB.

## Recommended Order

```matlab
addpath(genpath('analysis'));
cfg = default_config();

step01_preprocess_fmri(cfg);
step02_first_level_glm(cfg);
step03_define_fit_input_models(cfg);
step04_loo_input_model_validation(cfg);
step05_estimate_b_matrix_peb(cfg);
step06_compute_variance_explained(cfg);
step07_parameter_recovery_mc(cfg);
step08_detection_sensitivity_power(cfg);
```

Some steps require fitted outputs from earlier steps and can take a long time.

Preprocessing and first-level GLM steps rewrite large derivative files, so `run_pipeline` leaves them disabled unless `cfg.run_preprocessing` and `cfg.run_glm` are set to `true`.

```matlab
cfg.run_all = true;
cfg.run_preprocessing = true;
cfg.run_glm = true;
cfg.preprocessing.run_onset_extraction = false;
run_pipeline(cfg);
```

Set `cfg.preprocessing.run_onset_extraction = true` only when `correct_stimu_time.m` is available on the MATLAB path.

## Notes

- Paths are centralized in `config/default_config.m`.
- Large derived files are written to `cfg.output_dir`.
- Subject exclusions are configured in `cfg.exclude_subject_indices`.
- The DCM input-model space is the three-model space described in the 20260507 manuscript.
- `step05_estimate_b_matrix_peb` also writes B-matrix family posterior probabilities for the six family definitions used in the manuscript.
