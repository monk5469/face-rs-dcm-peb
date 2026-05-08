function cfg = default_config()
%DEFAULT_CONFIG Central configuration for public analysis scripts.

this_file = mfilename('fullpath');
analysis_dir = fileparts(fileparts(this_file));
repo_root = fileparts(analysis_dir);

cfg.repo_root = repo_root;
cfg.data_dir = fullfile(repo_root, 'data');
cfg.output_dir = fullfile(repo_root, 'outputs');

% The original local workflow used workspace_data.mat to define objects_dir
% and target_dir. Keep this configurable because those paths are not public.
cfg.workspace_file = fullfile(cfg.data_dir, 'workspace_data.mat');
cfg.bids_session_name = 'ses-mri';
cfg.functional_dir_name = 'func';
cfg.anatomical_dir_name = 'anat';
cfg.first_level_dir_name = 'Imm_Del_with6kernel';
cfg.dcm_dir_name = 'Imm_Del_DCM';
cfg.roi_dir_name = 'Imm_Del_ROI_adjust_eye(9)_Tcon_2_with6kernel_manual';

% The manuscript reports exclusion of two participants without both rOFA and
% rFFA VOIs. These are indices in the local objects_dir/target_dir arrays.
cfg.exclude_subject_indices = [7 12];

cfg.session_count = 1;
cfg.scale_modulatory_inputs = true;
cfg.region_names = {'rEVC', 'rOFA', 'rFFA', 'lEVC', 'lOFA', 'lFFA'};
cfg.voi_indices = [2 3]; % rOFA and rFFA

cfg.model_names = {'Model1_C1_All', 'Model2_C2_FaceScr', 'Model3_C3_NoRep'};
cfg.model_labels = {'M1 all visual', 'M2 face/scrambled', 'M3 initial only'};
cfg.selected_input_model = 2;

cfg.use_parfor = true;
cfg.parallel_workers = 10;
cfg.dcm_te = 0.03;
cfg.dcm_maxit = 128;

cfg.preprocessing.remove_dummy_scans = true;
cfg.preprocessing.run_realign = true;
cfg.preprocessing.run_slice_timing = true;
cfg.preprocessing.run_coregistration = true;
cfg.preprocessing.run_segmentation = true;
cfg.preprocessing.run_normalise_smooth = true;
cfg.preprocessing.run_onset_extraction = false;
cfg.preprocessing.dummy_scan_pattern = '*0000[12].nii';
cfg.preprocessing.functional_pattern = '^sub.*nii';
cfg.preprocessing.realign_pattern = '^rsub.*nii';
cfg.preprocessing.slice_timed_pattern = '^arsub.*nii';
cfg.preprocessing.mean_functional_pattern = '^meansub.*nii';
cfg.preprocessing.structural_pattern = '^sub.*mprage.*.nii';
cfg.preprocessing.tr = 2;
cfg.preprocessing.n_slices = 33;
cfg.preprocessing.ref_slice = 2;
cfg.preprocessing.slice_order = [1:2:33 2:2:32];
cfg.preprocessing.normalised_vox_functional = [3 3 3];
cfg.preprocessing.normalised_vox_structural = [1 1 1];
cfg.preprocessing.normalise_interp = 7;
cfg.preprocessing.normalise_wrap = [0 1 0];
cfg.preprocessing.smooth_fwhm = [6 6 6];

cfg.glm.condition_names = {'F_Init', 'F_Im', 'F_L', ...
    'U_Init', 'U_Im', 'U_L', 'S_Init', 'S_Im', 'S_L'};
cfg.glm.scans_pattern = '^s6war.*nii';
cfg.glm.motion_pattern = '^rp.*txt';
cfg.glm.onsets_file_name = 'trial_onsets.mat';
cfg.glm.run_first_level_spec = true;
cfg.glm.run_first_level_estimate = true;
cfg.glm.run_first_level_contrasts = true;
cfg.glm.run_second_level = false;
cfg.glm.second_level_dir_name = fullfile('GroupStats', 'Imm_Del');
cfg.glm.sessrep = 'none';
cfg.glm.hpf = 128;
cfg.glm.fmri_t = 16;
cfg.glm.fmri_t0 = 8;

cfg.peb_precision = 'single';
cfg.peb_fields_c = {'C'};
cfg.peb_fields_b = {'B'};
cfg.bma_probability_threshold = 0.95;

cfg.simulation.n_rep = 100;
cfg.simulation.snr_parameter_recovery = 6;
cfg.simulation.snr_detection = 2;
cfg.simulation.subject_jitter_sd = 0.10;
cfg.simulation.bad_subject_index = 12;
cfg.simulation.bad_subject_scale = 0.5;
cfg.simulation.target_connection_indices = [2 11 13 20];
cfg.simulation.rng_seed_parameter_recovery = 20251024;
cfg.simulation.rng_seed_detection = 20251030;
cfg.simulation.iterations = 1;
cfg.simulation.figure_flag = false;

cfg.run_all = false;
cfg.run_preprocessing = false;
cfg.run_glm = false;
end
