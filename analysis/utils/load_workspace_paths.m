function paths = load_workspace_paths(cfg)
%LOAD_WORKSPACE_PATHS Load and filter subject paths from workspace_data.mat.

if exist(cfg.workspace_file, 'file')
    S = load(cfg.workspace_file);
    required = {'objects_dir', 'target_dir'};
    for i = 1:numel(required)
        if ~isfield(S, required{i})
            error('Missing variable "%s" in %s.', required{i}, cfg.workspace_file);
        end
    end

    objects_dir = S.objects_dir;
    target_dir = S.target_dir;

    if isfield(S, 'sessions_dir')
        sessions_dir = S.sessions_dir;
    else
        sessions_dir = derive_sessions_dir(target_dir);
    end

    if isfield(S, 'target_structural_dir')
        target_structural_dir = S.target_structural_dir;
    else
        target_structural_dir = repmat({''}, size(target_dir));
    end
else
    [objects_dir, target_dir, target_structural_dir] = discover_workspace_from_data_dir(cfg);
    sessions_dir = derive_sessions_dir(target_dir);
end

exclude_idx = cfg.exclude_subject_indices;
exclude_idx = exclude_idx(exclude_idx >= 1 & exclude_idx <= numel(objects_dir));
objects_dir(exclude_idx) = [];
target_dir(exclude_idx) = [];
sessions_dir(exclude_idx) = [];
target_structural_dir(exclude_idx) = [];

paths.objects_dir = objects_dir;
paths.target_dir = target_dir;
paths.sessions_dir = sessions_dir;
paths.target_structural_dir = target_structural_dir;
paths.n_subjects = numel(target_dir);
end

function [objects_dir, target_dir, target_structural_dir] = discover_workspace_from_data_dir(cfg)
if ~exist(cfg.data_dir, 'dir')
    error('Neither workspace file nor data directory exists: %s', cfg.data_dir);
end

objects_dir = dir(cfg.data_dir);
objects_dir = objects_dir([objects_dir.isdir]);
objects_dir = objects_dir(~ismember({objects_dir.name}, {'.', '..'}));
objects_dir = objects_dir(~startsWith({objects_dir.name}, '.'));

target_dir = cell(1, numel(objects_dir));
target_structural_dir = cell(1, numel(objects_dir));
for s = 1:numel(objects_dir)
    subject_root = fullfile(cfg.data_dir, objects_dir(s).name);
    func_dir = fullfile(subject_root, cfg.bids_session_name, cfg.functional_dir_name);
    anat_dir = fullfile(subject_root, cfg.bids_session_name, cfg.anatomical_dir_name);

    if exist(func_dir, 'dir')
        target_dir{s} = func_dir;
    else
        target_dir{s} = subject_root;
    end

    if exist(anat_dir, 'dir')
        target_structural_dir{s} = anat_dir;
    else
        target_structural_dir{s} = fullfile(subject_root, cfg.anatomical_dir_name);
    end
end
end

function sessions_dir = derive_sessions_dir(target_dir)
sessions_dir = cell(size(target_dir));
for s = 1:numel(target_dir)
    entries = dir(target_dir{s});
    entries = entries([entries.isdir]);
    entries = entries(~ismember({entries.name}, {'.', '..', 'bold'}));
    sessions_dir{s} = entries;
end
end
