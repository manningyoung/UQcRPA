function [ unk ] = transform_bloch(datadir,unitary_matrix,wann_bands,kpoints)
% kpoints is only needed for length(kpoints)...

num_wann = length(unitary_matrix{1});
num_kpts = length(kpoints);

% Read in the reference Bloch wavefunctions at all k-points.
% The k-point index is added to the cell array as u0{ik,iband}(ir).
fprintf('Found %d files to read...\n',num_kpts);
for ik = 1:num_kpts
    fname = sprintf('%sUNK%05d.1',datadir,ik); % assuming all spin up (.1)
    fprintf('Reading %s\n',fname);
    temp = read_bloch(fname);
    u0(ik,:) = temp(wann_bands); % isolate Wannier bands
end

% Marzari & Vanderbilt (1997), Eq. (59)
fprintf('Transforming...\n',num_kpts);
num_grid = length(u0{1,1});
for ik = 1:num_kpts
    for n = 1:num_wann
        mysum = zeros(1,num_grid); % note: vectorized sum
        for m = 1:num_wann % sum over wanniers m
            mysum = mysum + unitary_matrix{ik}(m,n)*u0{ik,m};
        end
        unk{ik,n} = mysum;
    end
end
end          