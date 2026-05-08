function F = extract_free_energy_matrix(GCM)
%EXTRACT_FREE_ENERGY_MATRIX Extract subject x model free energies.

[n_sub, n_mod] = size(GCM);
F = nan(n_sub, n_mod);

for s = 1:n_sub
    for m = 1:n_mod
        if isempty(GCM{s, m})
            continue;
        end

        if ischar(GCM{s, m}) || isstring(GCM{s, m})
            D = load(GCM{s, m});
            if isfield(D, 'DCM')
                DCM = D.DCM;
            else
                continue;
            end
        else
            DCM = GCM{s, m};
        end

        if isfield(DCM, 'F')
            F(s, m) = DCM.F;
        end
    end
end
end

