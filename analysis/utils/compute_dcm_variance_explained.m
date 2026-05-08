function ve = compute_dcm_variance_explained(DCM)
%COMPUTE_DCM_VARIANCE_EXPLAINED Return VE by ROI for one fitted DCM.

if ~isfield(DCM, 'y') || ~isfield(DCM, 'R')
    ve = nan(1, DCM.n);
    return;
end

y_fit = DCM.y;
y_obs = DCM.y + DCM.R;
n_roi = size(y_obs, 2);
ve = nan(1, n_roi);

for r = 1:n_roi
    residual_ss = sum((y_obs(:, r) - y_fit(:, r)).^2);
    total_ss = sum((y_obs(:, r) - mean(y_obs(:, r))).^2);
    if total_ss > 0
        ve(r) = 1 - residual_ss / total_ss;
    end
end
end

