% UQcRPA - main job script
% See README for instructions
%
% TODO: clean up

datadir = '~/Downloads/data/';
wann_bands = 1:8;
iband = 1;
jband = 1;

%% CONSTRUCTION OF THE SCREENED INTERACTION W_{GG'}(q)

tic;

% Expected filenames
f_chi0mat = 'chi0mat.h5';
f_chimat = 'chimat.h5';
f_chi0mat_a = 'chi0mat.h5'; % FIXME: change back to chi0mat_a.h5
f_chimat_a = 'chimat.h5'; % FIXME: change back to chimat_a.h5
f_eps0mat = 'eps0mat.h5';
f_epsmat = 'epsmat.h5';

% Read polarisability (chi) matrices
chi0mat = h5read([datadir f_chi0mat],'/mats/matrix');
chimat = h5read([datadir f_chimat],'/mats/matrix');
chi0mat_a = h5read([datadir f_chi0mat_a],'/mats/matrix');
chimat_a = h5read([datadir f_chimat_a],'/mats/matrix');

chi0mat_a = zeros(size(chi0mat)); % FIXME: remove after testing
chimat_a = zeros(size(chimat)); % FIXME: remove after testing

% And their true dimensions, since the matrices are padded up to max(nmtx)
nmtx0 = h5read([datadir f_chi0mat],'/eps_header/gspace/nmtx');
nmtx = h5read([datadir f_chimat],'/eps_header/gspace/nmtx');
nmtx0_a = h5read([datadir f_chi0mat_a],'/eps_header/gspace/nmtx');
nmtx_a = h5read([datadir f_chimat_a],'/eps_header/gspace/nmtx');

if sum([nmtx0; nmtx] ~= [nmtx0_a; nmtx_a]) ~= 0
    error('P and P_a must be the same size.');
end

% Read Coulomb interaction from epsmat, since it's not initialised in chimat
vcoul0 = h5read([datadir f_eps0mat],'/eps_header/gspace/vcoul');
vcoul = h5read([datadir f_epsmat],'/eps_header/gspace/vcoul');

% Construct inverse dielectric matrix
[epsmat, vcoul_full] = calc_epsilon(chi0mat,chimat,chi0mat_a,chimat_a,nmtx0,nmtx,vcoul0,vcoul);

% Construct the screened interaction matrix
[W_matrix, Wp_diag] = calc_screened_matrix(epsmat,vcoul_full,[nmtx0; nmtx]);

%% CONSTUCTION OF THE TRANSFORMED BLOCH STATES u_{nk}(r)

% Read in the unitary transformation matrix U^(k)
[unitary_matrix, kpoints] = read_unitary([datadir 'silicon_u.mat']);

% Read in all UNKp.s files and apply the transformation
if exist([datadir 'unk.mat']) == 2
    fprintf('Loading precalculated unk from file...');
    load([datadir 'unk.mat']);
    fprintf('Done.\n');
else
    unk = transform_bloch(datadir,unitary_matrix,wann_bands,kpoints);
end

%% CALCULATION OF THE AUXILIARY FUNCTIONS F_{ij}(G,q)

% Read in the G-vectors corresponding to the matrix elements above
[gvecs, qpoints] = read_gvecs([datadir 'epsilon.log']);

% And the dimensions of the FFT
Nfft = double(h5read([datadir 'epsmat.h5'],'/mf_header/gspace/FFTgrid'));

% And the weights for the k-sum
kweights = read_kweights(datadir);

% Compute F_{ii}, F_{jj}
[Fii, gindex, ekin] = calc_aux(unk,kpoints,kweights,qpoints,Nfft,gvecs,iband,iband);
if iband ~= jband
    Fjj = calc_aux(unk,kpoints,kweights,qpoints,Nfft,gvecs,jband,jband);
else
    Fjj = Fii;
end
   
%% CALCULATION OF THE MATRIX ELEMENT U_{ij}

% W_matrix is ordered wrt the epsilon G-space while Fii, Fjj are not;
% we require gindexii and gindexjj to map to the correct index
Uij = calc_matrix_element(W_matrix,Fii,Fjj,gindex,[nmtx0; nmtx]);
fprintf('Calculation finished!\n');
Uij
toc

%% PLOTTING - MESSY
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
