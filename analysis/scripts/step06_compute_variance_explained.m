function step06_compute_variance_explained(cfg)
%STEP06_COMPUTE_VARIANCE_EXPLAINED Summarize DCM variance explained by ROI.

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
[n_sub, n_mod] = size(GCM);

rows = {};
for s = 1:n_sub
    for m = 1:n_mod
        DCM = GCM{s, m};
        ve = compute_dcm_variance_explained(DCM);
        for r = 1:numel(ve)
            rows(end + 1, :) = {s, input_models{m}.name, r, DCM.Y.name{r}, ve(r)}; %#ok<AGROW>
        end
    end
end

T = cell2table(rows, 'VariableNames', ...
    {'SubjectIndex', 'ModelName', 'ROIIndex', 'ROIName', 'VarianceExplained'});

out_dir = fullfile(cfg.output_dir, 'model_adequacy');
ensure_dir(out_dir);
writetable(T, fullfile(out_dir, 'variance_explained_by_roi.csv'));
save(fullfile(out_dir, 'variance_explained_by_roi.mat'), 'T', 'cfg');
end
