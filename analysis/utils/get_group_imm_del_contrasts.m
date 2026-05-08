function contrasts = get_group_imm_del_contrasts()
%GET_GROUP_IMM_DEL_CONTRASTS Return second-level Imm-Del contrasts.

contrasts.f_names = {'Imm_Del_repetition'};
contrasts.f_weights = { ...
    [0 1 -1 0 0 0 0 0 0; ...
     0 0 0 0 1 -1 0 0 0; ...
     0 0 0 0 0 0 0 1 -1]};

contrasts.t_names = {'main_effects_RS_faces'};
contrasts.t_weights = {[0 1 -1 0 1 -1 0 1 -1]};
end
