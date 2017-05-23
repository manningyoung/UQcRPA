% ==============================================================================
%
%   UQcRPA (0.1) main job script
%   See README for instructions
%
% ==============================================================================

%% User input

datadir = '~/Downloads/data/';
seedname = 'SrVO3';
wann_bands = 1:8;
iband = 1;
jband = 1;

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
[W_matrix, Wp_diag] = calc_screened_matrix(epsmat,vcoul_full,[nmtx0; nmtx]);

%% Calculation of the transformed Bloch states

% Read unitary transformation matrix
[unitary_matrix, kpoints] = read_unitary([datadir seedname '_u.mat']);

% Apply transformation to the reference Bloch states
if exist([datadir 'unk.mat']) == 2
    fprintf('Loading precalculated unk from file...');
    load([datadir 'unk.mat']);
    fprintf('Done.\n');
else
    unk = transform_bloch(datadir,unitary_matrix,wann_bands,kpoints);
end

%% Calculation of the auxiliary functions

% Read in the G-vectors corresponding to the screened matrix elements
[gvecs, qpoints] = read_gvecs([datadir 'epsilon.log']);

% And the dimensions of the FFT
Nfft = double(h5read([datadir 'epsmat.h5'],'/mf_header/gspace/FFTgrid'));

% And the k-point weights
kweights = read_kweights(datadir);

% Compute diagonal aux. functions F_{i,i} and F_{j,j} (if necessary)
[Fii, gindex, ekin] = calc_aux(unk,kpoints,kweights,qpoints,Nfft,gvecs,iband,iband);

if iband ~= jband
    Fjj = calc_aux(unk,kpoints,kweights,qpoints,Nfft,gvecs,jband,jband);
else
    Fjj = Fii;
end
   
%% Calculation of the matrix element

Uij = calc_matrix_element(W_matrix,Fii,Fjj,gindex,[nmtx0; nmtx]);
t_elapsed = toc;
fprintf('=========================\n');
fprintf('U_{%d,%d} = %.2f%+.2fi Ry\n',iband,jband,real(Uij),imag(Uij));
fprintf('=========================\n');
fprintf('Program completed successfully in %.1fs!\n',t_elapsed);

%% Misc. plotting/debugging

% close all
% 
% figure
% [sqrtekin, sortidx] = sort(sqrt(ekin{1}));
% plot(sqrtekin,abs(Fii{1}(sortidx)).^2,':');
% xlabel('|q+G|')
% ylabel('|F_{11}(G,q=0)|^2')
% xlim([0 10])
% hold on
% plot(sqrtekin,abs(Fii{1}(sortidx)).^2,'x');
%

% Plot the diagonal of the polarization potenal W^p=(W-v)(q=0)
figure
for ig=2:length(gvecs{1})
    mynorm(ig-1) = norm(qpoints(1,:)+gvecs{1}(ig,:)); % |q+G|
end
[sortmynorm,sortidx] = sort(mynorm); % |q+G| sorted
plot(sortmynorm,abs(Wp_diag{1}(sortidx)),'x'); % plot vs W^p with same sorting
xlabel('|q+G|')
ylabel('|W^p_{GG}(q=0)|')
grid on
hold on
plot(sortmynorm,abs(Wp_diag{1}(sortidx)),':');
