function step02_first_level_glm(cfg)
%STEP02_FIRST_LEVEL_GLM Run Imm-Del first-level and optional group GLM.

if nargin < 1 || isempty(cfg)
    cfg = default_config();
end

paths = load_workspace_paths(cfg);
spm('defaults', 'FMRI');

if cfg.glm.run_first_level_spec
    specify_first_level_models(paths, cfg);
end

if cfg.glm.run_first_level_estimate
    estimate_first_level_models(paths, cfg);
end

if cfg.glm.run_first_level_contrasts
    define_first_level_contrasts(paths, cfg);
end

if cfg.glm.run_second_level
    run_second_level_anova(paths, cfg);
end
end

function specify_first_level_models(paths, cfg)
for s = 1:paths.n_subjects
    session_dirs = paths.sessions_dir{s};
    out_dir = fullfile(paths.target_dir{s}, 'bold', cfg.first_level_dir_name);
    ensure_dir(out_dir);

    onset_data = load(fullfile(paths.target_dir{s}, 'bold', cfg.glm.onsets_file_name), 'trial_onsets');
    trial_onsets = onset_data.trial_onsets;

    clear matlabbatch
    matlabbatch{1}.spm.stats.fmri_spec.dir = {out_dir};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = cfg.preprocessing.tr;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = cfg.glm.fmri_t;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = cfg.glm.fmri_t0;

    [scan_list, motion_file, onsets, nscan] = concatenate_session_inputs( ...
        paths.target_dir{s}, session_dirs, trial_onsets, cfg);

    matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = scan_list;
    for c = 1:numel(cfg.glm.condition_names)
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(c).name = cfg.glm.condition_names{c};
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(c).onset = onsets{c};
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(c).duration = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(c).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(c).pmod = struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(c).orth = 1;
    end

    matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = {motion_file};
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = cfg.glm.hpf;
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

    spm_jobman('run', matlabbatch);
    spm_fmri_concatenate(fullfile(out_dir, 'SPM.mat'), nscan);
end
end

function [scan_list, motion_file, onsets, nscan] = concatenate_session_inputs(target_dir, session_dirs, trial_onsets, cfg)
scan_list = {};
R = [];
onsets = cell(1, numel(cfg.glm.condition_names));
for c = 1:numel(onsets)
    onsets{c} = [];
end

time_so_far = 0;
nscan = zeros(1, numel(session_dirs));
for i = 1:numel(session_dirs)
    session_path = fullfile(target_dir, session_dirs(i).name);
    scans = cellstr(spm_select('FPList', session_path, cfg.glm.scans_pattern));
    scan_list = [scan_list; scans]; %#ok<AGROW>

    rp_file = spm_select('FPList', session_path, cfg.glm.motion_pattern);
    R = [R; load(rp_file)]; %#ok<AGROW>

    for c = 1:numel(cfg.glm.condition_names)
        onsets{c} = [onsets{c}; trial_onsets{1, i}{1, c} + time_so_far];
    end

    nscan(i) = numel(scans);
    time_so_far = time_so_far + nscan(i) * cfg.preprocessing.tr;
end

motion_file = fullfile(target_dir, 'bold', cfg.first_level_dir_name, 'rp_all_sess.mat');
save(motion_file, 'R');
end

function estimate_first_level_models(paths, cfg)
for s = 1:paths.n_subjects
    spm_file = fullfile(paths.target_dir{s}, 'bold', cfg.first_level_dir_name, 'SPM.mat');
    S = load(spm_file, 'SPM');
    spm_spm(S.SPM);
end
end

function define_first_level_contrasts(paths, cfg)
contrasts = get_imm_del_contrasts();

for s = 1:paths.n_subjects
    clear matlabbatch
    spm_file = spm_select('FPList', fullfile(paths.target_dir{s}, 'bold', cfg.first_level_dir_name), '^SPM.mat');
    matlabbatch{1}.spm.stats.con.spmmat = cellstr(spm_file);

    for i = 1:numel(contrasts.f_names)
        matlabbatch{1}.spm.stats.con.consess{i}.fcon.name = contrasts.f_names{i};
        matlabbatch{1}.spm.stats.con.consess{i}.fcon.convec = contrasts.f_weights{i};
        matlabbatch{1}.spm.stats.con.consess{i}.fcon.sessrep = cfg.glm.sessrep;
    end

    offset = numel(contrasts.f_names);
    for i = 1:numel(contrasts.t_names)
        matlabbatch{1}.spm.stats.con.consess{i + offset}.tcon.name = contrasts.t_names{i};
        matlabbatch{1}.spm.stats.con.consess{i + offset}.tcon.convec = contrasts.t_weights{i};
        matlabbatch{1}.spm.stats.con.consess{i + offset}.tcon.sessrep = cfg.glm.sessrep;
    end

    matlabbatch{1}.spm.stats.con.delete = 1;
    spm_jobman('run', matlabbatch);
end
end

function run_second_level_anova(paths, cfg)
output_dir = fullfile(cfg.data_dir, cfg.glm.second_level_dir_name);
ensure_dir(output_dir);

clear matlabbatch
matlabbatch{1}.spm.stats.factorial_design.dir = {output_dir};
for s = 1:paths.n_subjects
    contrasts = cellstr(spm_select('FPList', ...
        fullfile(paths.target_dir{s}, 'bold', cfg.first_level_dir_name), '^con.*nii'));
    matlabbatch{1}.spm.stats.factorial_design.des.anovaw.fsubject(s).scans = contrasts(1:9, 1);
    matlabbatch{1}.spm.stats.factorial_design.des.anovaw.fsubject(s).conds = 1:9;
end

matlabbatch{1}.spm.stats.factorial_design.des.anovaw.dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.anovaw.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.anovaw.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.anovaw.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
spm_jobman('run', matlabbatch);

S = load(fullfile(output_dir, 'SPM.mat'), 'SPM');
spm_spm(S.SPM);
define_second_level_contrasts(output_dir, cfg);
end

function define_second_level_contrasts(output_dir, cfg)
contrasts = get_group_imm_del_contrasts();

clear matlabbatch
matlabbatch{1}.spm.stats.con.spmmat = cellstr(spm_select('FPList', output_dir, '^SPM.*mat'));
for i = 1:numel(contrasts.f_names)
    matlabbatch{1}.spm.stats.con.consess{i}.fcon.name = contrasts.f_names{i};
    matlabbatch{1}.spm.stats.con.consess{i}.fcon.convec = contrasts.f_weights{i};
    matlabbatch{1}.spm.stats.con.consess{i}.fcon.sessrep = cfg.glm.sessrep;
end

offset = numel(contrasts.f_names);
for i = 1:numel(contrasts.t_names)
    matlabbatch{1}.spm.stats.con.consess{i + offset}.tcon.name = contrasts.t_names{i};
    matlabbatch{1}.spm.stats.con.consess{i + offset}.tcon.convec = contrasts.t_weights{i};
    matlabbatch{1}.spm.stats.con.consess{i + offset}.tcon.sessrep = cfg.glm.sessrep;
end

matlabbatch{1}.spm.stats.con.delete = 1;
spm_jobman('run', matlabbatch);
end
