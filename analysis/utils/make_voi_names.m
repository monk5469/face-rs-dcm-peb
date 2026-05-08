function voi_names = make_voi_names(cfg)
%MAKE_VOI_NAMES Build SPM VOI file stems for all regions and sessions.

voi_names = cell(numel(cfg.region_names), cfg.session_count);
for v = 1:numel(cfg.region_names)
    for ses = 1:cfg.session_count
        voi_names{v, ses} = sprintf('VOI_%s_%d', cfg.region_names{v}, ses);
    end
end
end

