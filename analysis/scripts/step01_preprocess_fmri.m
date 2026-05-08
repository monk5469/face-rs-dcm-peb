function step01_preprocess_fmri(cfg)
%STEP01_PREPROCESS_FMRI Run SPM preprocessing for the fMRI data.

if nargin < 1 || isempty(cfg)
    cfg = default_config();
end

paths = load_workspace_paths(cfg);
spm('defaults', 'FMRI');

for s = 1:paths.n_subjects
    session_dirs = paths.sessions_dir{s};

    if cfg.preprocessing.remove_dummy_scans
        move_dummy_scans(paths.target_dir{s}, session_dirs, cfg);
    end

    if cfg.preprocessing.run_realign
        run_realign(paths.target_dir{s}, session_dirs, cfg);
    end

    if cfg.preprocessing.run_slice_timing
        run_slice_timing(paths.target_dir{s}, session_dirs, cfg);
    end

    if cfg.preprocessing.run_coregistration
        run_coregistration(paths.target_dir{s}, paths.target_structural_dir{s}, session_dirs, cfg);
    end

    if cfg.preprocessing.run_segmentation || cfg.preprocessing.run_normalise_smooth
        sn = run_segmentation(paths.target_structural_dir{s}, cfg);
    else
        sn = [];
    end

    if cfg.preprocessing.run_normalise_smooth
        run_normalise_smooth(paths.target_dir{s}, session_dirs, sn, cfg);
    end

    if cfg.preprocessing.run_onset_extraction
        run_onset_extraction(paths, s, cfg);
    end
end
end

function move_dummy_scans(target_dir, session_dirs, cfg)
unused_dir = fullfile(target_dir, 'bold', 'unused_scans');
ensure_dir(unused_dir);

for i = 1:numel(session_dirs)
    session_path = fullfile(target_dir, session_dirs(i).name);
    scans = dir(fullfile(session_path, cfg.preprocessing.dummy_scan_pattern));
    for j = 1:numel(scans)
        movefile(fullfile(session_path, scans(j).name), unused_dir);
    end
end
end

function run_realign(target_dir, session_dirs, cfg)
P = cell(numel(session_dirs), 1);
for i = 1:numel(session_dirs)
    P{i} = spm_select('FPList', fullfile(target_dir, session_dirs(i).name), ...
        cfg.preprocessing.functional_pattern);
end

flags = struct('rtn', 1);
spm_realign(P, flags);
spm_reslice(P, flags);
end

function run_slice_timing(target_dir, session_dirs, cfg)
PF = cell(numel(session_dirs), 1);
for i = 1:numel(session_dirs)
    PF{i} = spm_select('FPList', fullfile(target_dir, session_dirs(i).name), ...
        cfg.preprocessing.realign_pattern);
end

timing = [cfg.preprocessing.tr / cfg.preprocessing.n_slices, ...
    cfg.preprocessing.tr / cfg.preprocessing.n_slices];
spm_slice_timing(PF, cfg.preprocessing.slice_order, ...
    cfg.preprocessing.ref_slice, timing);
end

function run_coregistration(target_dir, structural_dir, session_dirs, cfg)
PS = spm_select('FPList', structural_dir, cfg.preprocessing.structural_pattern);
PF = spm_select('FPList', fullfile(target_dir, session_dirs(1).name), ...
    cfg.preprocessing.mean_functional_pattern);

flags = [];
out = spm_coreg(PS, PF, flags);
M = spm_matrix(out);
MM = spm_get_space(PS);
spm_get_space(PS, M \ MM);
end

function sn = run_segmentation(structural_dir, cfg)
PS = spm_select('FPList', structural_dir, cfg.preprocessing.structural_pattern);
sn_file = sprintf('%s_seg_sn.mat', spm_str_manip(PS, 'sd'));

if exist(sn_file, 'file')
    S = load(sn_file, 'sn');
    sn = S.sn;
    return;
end

out = spm_preproc(PS);
[sn, isn] = spm_prep2sn(out);
save(sn_file, 'sn');
save(sprintf('%s_seg_inv_sn.mat', spm_str_manip(PS, 'sd')), 'isn');

opts.biascor = 1;
opts.GM = [1 0 1];
opts.WM = [1 0 0];
opts.CSF = [1 0 0];
opts.cleanup = 1;
spm_preproc_write(sn, opts);
end

function run_normalise_smooth(target_dir, session_dirs, sn, cfg)
flags.interp = cfg.preprocessing.normalise_interp;
flags.wrap = cfg.preprocessing.normalise_wrap;
flags.vox = cfg.preprocessing.normalised_vox_functional;

if cfg.use_parfor
    start_parallel_pool(cfg.parallel_workers);
end

parfor sess = 1:numel(session_dirs)
    PF = spm_select('FPList', fullfile(target_dir, session_dirs(sess).name), ...
        cfg.preprocessing.slice_timed_pattern);
    for m = 1:size(PF, 1)
        VO = spm_write_sn(PF(m, :), sn, flags);
        [pth, nam, ext] = fileparts(PF(m, :));
        spm_smooth(VO, fullfile(pth, ['s6w' nam ext]), cfg.preprocessing.smooth_fwhm);
    end
end
end

function run_onset_extraction(paths, subject_index, cfg)
if exist('correct_stimu_time', 'file') ~= 2
    error(['cfg.preprocessing.run_onset_extraction is true, but ' ...
        'correct_stimu_time.m is not on the MATLAB path.']);
end

target_dir = paths.target_dir{subject_index};
session_dirs = paths.sessions_dir{subject_index};
fname = dir(fullfile(target_dir, '*.tsv'));
trial_onsets = correct_stimu_time(cfg.data_dir, paths.objects_dir(subject_index).name, ...
    target_dir, session_dirs, fname, subject_index);
save(fullfile(target_dir, 'bold', cfg.glm.onsets_file_name), 'trial_onsets');
end
