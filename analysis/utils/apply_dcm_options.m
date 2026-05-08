function DCM = apply_dcm_options(DCM, cfg)
%APPLY_DCM_OPTIONS Apply the fMRI DCM options used throughout the analysis.

DCM.options.nonlinear = 0;
DCM.options.two_state = 0;
DCM.options.stochastic = 0;
DCM.options.centre = 1;
DCM.options.nograph = 1;
DCM.options.maxnodes = 8;
DCM.options.maxit = cfg.dcm_maxit;
DCM.options.hidden = [];
DCM.options.induced = 0;
DCM.M.options = struct();
end

