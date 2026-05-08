# Face Repetition Suppression DCM-PEB Analysis Code

This repository contains MATLAB/SPM12 analysis scripts for the fMRI DCM-PEB analyses reported in:

Face Repetition Suppression Reflects Condition-Dependent Reconfiguration of Occipitotemporal Effective Connectivity.

The standardized analysis workflow lives in `analysis/`.

## Purpose and Attribution

This repository supports academic transparency, peer review, and reproducible analysis for the accompanying face repetition suppression DCM-PEB study.

Copyright (c) 2026 Yunnan Minzu University, Jia Tang and Haifeng Wu.

版权所有 (c) 2026 云南民族大学 唐嘉、吴海锋。

Parts of the analysis workflow and code structure were adapted and modified from analysis code associated with Lee et al. (2022), which examined effects of face repetition on ventral visual stream connectivity using DCM of fMRI data. The present repository extends that framework for the current input-model comparison, condition-specific B-matrix PEB/BMA, family-level B-matrix posterior probabilities, and simulation-based validation analyses.

## Requirements

- MATLAB with the Statistics and Machine Learning Toolbox
- SPM12 on the MATLAB path
- Parallel Computing Toolbox for Monte Carlo simulations
- OpenNeuro dataset `ds000117`, fMRI component
- Preprocessed first-level SPM directories and VOI files generated from the ROI-localization pipeline described in the manuscript

## Quick Start

1. Edit `analysis/config/default_config.m` so that paths match your local SPM first-level and VOI derivatives. If `workspace_data.mat` is unavailable, the loader falls back to scanning `cfg.data_dir` using the configured BIDS-style session names.
2. From MATLAB, run:

```matlab
addpath(genpath('analysis'));
cfg = default_config();
run_pipeline(cfg);
```

By default, `run_pipeline` executes only the DCM input-model definition and fitting step, because later steps depend on fitted outputs and can be computationally expensive. Individual steps can be run directly from `analysis/scripts/`.

## Repository Layout

- `analysis/config/default_config.m`: central configuration for data paths, subject exclusions, model names, and simulation settings.
- `analysis/scripts/step01_preprocess_fmri.m`: SPM preprocessing pipeline, including dummy-scan removal, realignment, slice timing, coregistration, segmentation, normalization, and smoothing.
- `analysis/scripts/step02_first_level_glm.m`: first-level Imm-Del GLM specification, concatenation, estimation, contrasts, and optional second-level ANOVA.
- `analysis/scripts/step03_define_fit_input_models.m`: defines and fits the three candidate DCM input architectures.
- `analysis/scripts/step04_loo_input_model_validation.m`: leave-one-out validation for input-model comparison.
- `analysis/scripts/step05_estimate_b_matrix_peb.m`: group-level PEB/BMR/BMA for condition-specific B-matrix modulation and family posterior probabilities.
- `analysis/scripts/step06_compute_variance_explained.m`: variance-explained model adequacy summary.
- `analysis/scripts/step07_parameter_recovery_mc.m`: Monte Carlo parameter-recovery simulation.
- `analysis/scripts/step08_detection_sensitivity_power.m`: simulation-based detection-sensitivity analysis.
- `analysis/utils/`: helper functions shared by the scripts.
- `CODE_COVERAGE_FOR_MANUSCRIPT.md`: mapping between manuscript analyses and repository scripts.
- `CODE_METADATA.md`: code purpose, authorship, dependency, and licensing notes.
- `LICENSE`: MIT open-source license.
- `NOTICE.md`: copyright, purpose, and attribution notice.

## Data Availability

Raw data are available from OpenNeuro dataset `ds000117`. Large derived files such as fitted GCM/DCM `.mat` files, figures, and simulation outputs are not included in this code repository.

## License

This code is released under the MIT License to support open, reusable, and reproducible research. See `LICENSE`, `NOTICE.md`, and `CODE_METADATA.md`.
