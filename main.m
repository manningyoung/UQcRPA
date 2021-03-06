function [ Uij ] = main(datadir,seedname,wann_bands,iband,jband)
% ==============================================================================
%
%   UQcRPA main program
%   See README for details
%
% ==============================================================================

if nargin < 5
    error('Not enough input arguments.');
end

%% Construction of the screened interaction matrix

tic;

% Read polarisability (chi) matrices
chi0mat = h5read([datadir 'chi0mat.h5'],'/mats/matrix');
chimat = h5read([datadir 'chimat.h5'],'/mats/matrix');
chi0mat_a = h5read([datadir 'chi0mat.h5'],'/mats/matrix'); % FIXME: change back to chi0mat_a.h5
chimat_a = h5read([datadir 'chimat.h5'],'/mats/matrix'); % FIXME: change back to chimat_a.h5

chi0mat_a = zeros(size(chi0mat)); % FIXME: remove after testing
chimat_a = zeros(size(chimat)); % FIXME: remove after testing

% Read matrix dimensions, since they're padded up to max(nmtx)
nmtx0 = h5read([datadir 'chi0mat.h5'],'/eps_header/gspace/nmtx');
nmtx = h5read([datadir 'chimat.h5'],'/eps_header/gspace/nmtx');
nmtx0_a = h5read([datadir 'chi0mat.h5'],'/eps_header/gspace/nmtx'); % FIXME: change back to chi0mat_a.h5
nmtx_a = h5read([datadir 'chimat.h5'],'/eps_header/gspace/nmtx'); % FIXME: change back to chimat_a.h5

if sum([nmtx0; nmtx] ~= [nmtx0_a; nmtx_a]) ~= 0
    error('P and P_a must be the same size.');
end

% Read Coulomb interaction from epsmat, since it's not initialised in chimat
vcoul0 = h5read([datadir 'eps0mat.h5'],'/eps_header/gspace/vcoul');
vcoul = h5read([datadir 'epsmat.h5'],'/eps_header/gspace/vcoul');

% Construct inverse dielectric matrix and screened interaction
[epsmat, vcoul_full] = calc_epsilon(chi0mat,chimat,chi0mat_a,chimat_a,nmtx0,nmtx,vcoul0,vcoul);
W_matrix = calc_screened_matrix(epsmat,vcoul_full,[nmtx0; nmtx]);

%% Calculation of the transformed Bloch states

% Read unitary transformation matrix
% Unnecessary if unk.mat exists below, but we still need kpoints...
[unitary_matrix, kpoints] = read_unitary([datadir seedname '_u.mat']);

% Apply transformation to the reference Bloch states
if exist([datadir 'unk.mat']) == 2
    fprintf('Loading precalculated unk from file...');
    unk_mat = load([datadir 'unk.mat']);
    unk = unk_mat.unk;
    fprintf('Done.\n');
else
    unk = transform_bloch(datadir,unitary_matrix,wann_bands,kpoints);
    fprintf('Saving unk to file for next time...');
    save([datadir 'unk.mat'],'unk');
    fprintf('Done.\n');
end

%% Calculation of the auxiliary functions

% Read in the G-vectors corresponding to the screened matrix elements
[gvecs, qpoints] = read_gvecs([datadir 'epsilon.log']);

% And the dimensions of the FFT
Nfft = double(h5read([datadir 'epsmat.h5'],'/mf_header/gspace/FFTgrid'));

% Assume the Wannier k-grid isn't symmetry reduced
kweights = ones(1,length(kpoints));

% Compute diagonal aux. functions F_{i,i} and F_{j,j} as necessary
Fii_str = sprintf('F_%d_%d.mat',iband,iband);
Fjj_str = sprintf('F_%d_%d.mat',jband,jband);

% Each F_iband_jband.mat file contains variables F and gindex
% F can then be loaded as Fii or Fjj as required
if exist([datadir Fii_str]) == 2
    fprintf('Loading precalculated Fii from file...');
    F_mat = load([datadir Fii_str]);
    Fii = F_mat.F;
    gindex = F_mat.gindex;
    fprintf('Done.\n');
else 
    [Fii, gindex] = calc_aux(unk,kpoints,kweights,qpoints,Nfft,gvecs,iband,iband);
    F = Fii;
    fprintf('Saving Fii to file for next time...');
    save([datadir Fii_str],'F','gindex');
    fprintf('Done.\n');
end

if exist([datadir Fjj_str]) == 2
    fprintf('Loading precalculated Fjj from file...');
    F_mat = load([datadir Fjj_str]);
    Fjj = F_mat.F;
    fprintf('Done.\n');
else 
    [Fjj, gindex] = calc_aux(unk,kpoints,kweights,qpoints,Nfft,gvecs,jband,jband);
    fprintf('Saving Fjj to file for next time...');
    F = Fjj;
    save([datadir Fjj_str],'F','gindex');
    fprintf('Done.\n');
end
   
%% Calculation of the matrix element

Uij = calc_matrix_element(W_matrix,Fii,Fjj,gindex,[nmtx0; nmtx]);
t_elapsed = toc;
fprintf('==========================\n');
fprintf('U_{%d,%d} = %.3f%+.3fi Ry\n',iband,jband,real(Uij),imag(Uij));
fprintf('==========================\n');
fprintf('Program completed successfully in %.1fs!\n',t_elapsed);

end
