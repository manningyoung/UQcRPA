function [ unk ] = transform_bloch(datadir,unitary_matrix,wann_bands,kpoints)
% Applies the unitary transformation matrix to the reference Bloch states in
% datadir for construction of the Wannier states. Only transforms the bands
% indexed by wann_bands (active subspace).
%
% FIXME: assumes non spin-polarised states (UNKp.s, s=1).

fprintf('Transforming Bloch states...\n');

num_wann = length(unitary_matrix{1});

% Read in the reference Bloch wavefunctions at all k-points.
% The k-point index is added to the cell array as u0{ik,iband}(ir).
fprintf('Found %d files to read.\n',length(kpoints));
for ik = 1:length(kpoints)
    fname = sprintf('%sUNK%05d.1',datadir,ik);
    temp = read_bloch(fname);
    u0(ik,:) = temp(wann_bands); % isolate Wannier bands
end

% Eq. (59), N. Marzari and D. Vanderbilt, Phys. Rev. B 56(20), 12847 (1997)
fprintf('Applying transformation...');
num_grid = length(u0{1,1});
for ik = 1:length(kpoints)
    for n = 1:num_wann
        mysum = zeros(1,num_grid); % note: vectorised sum
        for m = 1:num_wann % sum over wanniers m
            mysum = mysum + unitary_matrix{ik}(m,n)*u0{ik,m};
        end
        unk{ik,n} = mysum;
    end
end

fprintf('Done.\n');

end          