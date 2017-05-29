function [ Uij ] = calc_matrix_element(W_matrix,Fii,Fjj,gindex,nmtx)
% Calculates the screened Coulomb matrix element U_{i,j}.
%
% Fii and Fjj are calculated for all G, but the sum is performed over the subset
% in gindex due to the dielectric energy cutoff on W. W_matrix is ordered wrt
% the epsilon G-space (|q+G|^2) while Fii, Fjj are not; we require gindex to
% match the correct elements.

fprintf('Taking matrix element in the Wannier basis...');

Uij = 0;
for iq = 1:length(nmtx)
    for ig = 2:nmtx(iq) % FIXME: only body terms
        for igp = 2:nmtx(iq)
            if iq==34
                igp
            end
            Uij = Uij + conj(Fii{iq}(gindex{iq}(ig)))*W_matrix{iq}(ig,igp)*Fjj{iq}(gindex{iq}(igp));
        end
    end
end

fprintf('Done.\n');

end