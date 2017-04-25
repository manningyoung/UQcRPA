function [ W_matrix , Wp_diag ] = calc_screened_matrix(epsmat,vcoul,nmtx)
% Calculates the screened Coulomb interaction W via:
%   W_{GG'}(q) = epsilon_{GG'}^{-1}(q) * v(q+G').
% Assumes epsinv_matrix and vcoul are ordered wrt the epsilon G-space |q+G|^2.
% Also returns the diagonal elements of the polarization potential W^p = W-v.

% Deslippe et al. (2012), Eq. (12)
for iq = 1:length(nmtx)
    for ig = 2:nmtx(iq) % BODY ONLY
        for igp = 2:nmtx(iq)
            W_matrix{iq}(ig,igp) = epsmat{iq}(ig,igp)*vcoul{iq}(igp);
            if ig == igp
                Wp_diag{iq}(ig) = W_matrix{iq}(ig,ig)-vcoul{iq}(ig);
            end
        end
    end
end

end