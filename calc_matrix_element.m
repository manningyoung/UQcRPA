function [ Uij ] = calc_matrix_element(W_matrix,Fii,Fjj,gindex,nmtx)
% Calculates the matrix element U_{ij}.
% Fii, Fjj are calculated for all G but the sum below is only done over
% gvecs due to the cutoff on W.

Uij = 0;
for iq = 1:length(nmtx)
    for ig = 2:nmtx(iq) % BODY ONLY
        for igp = 2:nmtx(iq)
            Uij = Uij + conj(Fii{iq}(gindex{iq}(ig)))*W_matrix{iq}(ig,igp)*Fjj{iq}(gindex{iq}(igp));
        end
    end
end

end