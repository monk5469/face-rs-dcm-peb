function DCM = build_dcm_from_spm(SPM, xY, model_spec, cfg, session_index)
%BUILD_DCM_FROM_SPM Build one subject/session fMRI DCM from SPM and VOIs.

n_voi = numel(xY);
n_inputs = size(model_spec.contrast, 1);
driver_rows = model_spec.driver_rows;

DCM = struct();
DCM.xY = xY;
DCM.n = numel(DCM.xY);
DCM.v = numel(DCM.xY(1).u);
DCM.Y.dt = SPM.xY.RT;
DCM.Y.X0 = DCM.xY(1).X0;

for i = 1:DCM.n
    DCM.Y.y(:, i) = DCM.xY(i).u;
    DCM.Y.name{i} = DCM.xY(i).name;
end

DCM.Y.Q = spm_Ce(ones(1, DCM.n) * DCM.v);
DCM.U.dt = SPM.Sess(session_index).U(1).dt;

for ui = 1:n_inputs
    active_conditions = find(model_spec.contrast(ui, :) ~= 0);
    u = sparse(length(SPM.Sess(session_index).U(1).u) - 32, 1);

    for c = 1:numel(active_conditions)
        condition_index = active_conditions(c);
        weight = model_spec.contrast(ui, condition_index);
        u = u + weight * SPM.Sess(session_index).U(condition_index).u(33:end);
    end

    if cfg.scale_modulatory_inputs && ~ismember(ui, driver_rows)
        u_sd = std(full(u));
        if u_sd > 0
            u = u / u_sd;
        end
    end

    DCM.U.u(:, ui) = u;
    DCM.U.name{ui} = make_input_name(SPM, session_index, active_conditions, model_spec.contrast(ui, :));
end

DCM.delays = repmat(SPM.xY.RT, n_voi, 1) / 2;
DCM.TE = cfg.dcm_te;
DCM = apply_dcm_options(DCM, cfg);

DCM.a = ones(n_voi);
DCM.b = zeros(n_voi, n_voi, n_inputs);
DCM.c = zeros(n_voi, n_inputs);
DCM.d = zeros(DCM.n, DCM.n, 0);

for k = driver_rows
    DCM.c(:, k) = 1;
end

for k = 1:n_inputs
    if ~ismember(k, driver_rows)
        DCM.b(:, :, k) = ones(n_voi);
    end
end
end

function input_name = make_input_name(SPM, session_index, active_conditions, weights)
input_name = SPM.Sess(session_index).U(active_conditions(1)).name{1};
for c = 2:numel(active_conditions)
    condition_index = active_conditions(c);
    if weights(condition_index) == 1
        symbol = '+';
    else
        symbol = '-';
    end
    input_name = [input_name symbol SPM.Sess(session_index).U(condition_index).name{1}]; %#ok<AGROW>
end
end

