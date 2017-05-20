function [ W_matrix , Wp_diag ] = calc_screened_matrix(epsmat,vcoul,nmtx)
% Calculates the screened Coulomb interaction matrix according to
%
%   W_{GG'}(q) = \epsilon_{GG'}^{-1}(q)*v(q+G')
%
% Also returns the diagonal elements of the polarization potential W^p = W-v.

fprintf('Constructing screened interaction matrix...');

for iq = 1:length(nmtx)
    for ig = 2:nmtx(iq) % FIXME: only body terms
        for igp = 2:nmtx(iq)
            W_matrix{iq}(ig,igp) = epsmat{iq}(ig,igp)*vcoul{iq}(igp);
            if ig == igp
                Wp_diag{iq}(ig) = W_matrix{iq}(ig,ig)-vcoul{iq}(ig);
            end
        end
    end
end

fprintf('Done!\n');

end