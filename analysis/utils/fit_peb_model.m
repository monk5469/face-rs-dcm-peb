function [GCM_fit, PEB] = fit_peb_model(GCM, M, fields, cfg)
%FIT_PEB_MODEL Fit a PEB model using the best available local function.

if nargin < 4
    cfg = struct();
end
if ~isfield(cfg, 'use_parfor')
    cfg.use_parfor = false;
end
if ~isfield(cfg, 'simulation') || ~isfield(cfg.simulation, 'iterations')
    cfg.simulation.iterations = 1;
end

if exist('region_spm_dcm_peb_fit', 'file') == 2
    [GCM_fit, PEB] = region_spm_dcm_peb_fit( ...
        GCM, M, fields, cfg.use_parfor, cfg.simulation.iterations);
elseif exist('spm_dcm_peb_fit', 'file') == 2
    [GCM_fit, PEB] = spm_dcm_peb_fit(GCM, M, fields);
else
    GCM_fit = GCM;
    PEB = spm_dcm_peb(GCM_fit, M, fields);
end
end

