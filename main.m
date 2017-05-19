% UQcRPA
% Copyright (c) 2017, The University of Queensland
%
% Version 0.1 (May 2017)
% M. C. Young, A. C. Jacko, B. J. Powell
% MATLAB R2017a / BerkeleyGW 1.2 / wannier90 2.1.0
%
% Required inputs:
%   chimat.h5       [BerkeleyGW]
%   chi0mat.h5      [BerkeleyGW]
%   chimat_a.h5     [BerkeleyGW]
%   chi0mat_a.h5    [BerkeleyGW]
%   epsmat.h5       [BerkeleyGW]
%   eps0mat.h5      [BerkeleyGW]
%   epsilon.log     [BerkeleyGW]
%   seedname_u.mat  [wannier90]
%   {UNKXXXXX.1}    [wannier90]

% Set the location of the input directory (blank if PWD)
datadir = '~/Downloads/data/';

% Define the Wannier bands
wann_bands = 1:8;

% Define the matrix element U_{ij}
iband = 1;
jband = 1;

%% CONSTRUCTION OF THE SCREENED INTERACTION W_{GG'}(q)

% 
chimat = h5read([datadir 'chimat.h5'],'/mats/matrix');
chi0mat = h5read([datadir 'chimat.h5'],'/mats/matrix');
chimat_a = h5read([datadir 'chimat_a.h5'],'/mats/matrix');
chi0mat_a = h5read([datadir 'chimat_a.h5'],'/mats/matrix');
nmtx = 
nmtx0 = 
vcoul = h5read([datadir 'epsmat.h5'],'/eps_header/gspace/vcoul');
vcoul0 = h5read([datadir 'eps0mat.h5'],'/eps_header/gspace/vcoul');



epsmat = calc_epsilon(chi0mat,chimat,chi0mat_a,chimat_a,vcoul0,vcoul,nmtx0,nmtx);

% Construct the full W matrix
[W_matrix, Wp_diag] = calc_screened_matrix(epsmat_full,vcoul_full,nmtx_full);

%% CONSTUCTION OF THE TRANSFORMED BLOCH STATES u_{nk}(r)

% Read in the unitary transformation matrix U^(k)
[unitary_matrix, kpoints] = read_unitary([datadir 'silicon_u.mat']);

% Read in all UNKp.s files and apply the transformation
% unk = transform_bloch(datadir,unitary_matrix,wann_bands,kpoints);
% OR save workspace and load the precalculated u_{nk}:
load([datadir 'unk.mat']);

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
Uij = calc_matrix_element(W_matrix,Fii,Fjj,gindex,nmtx_full);
fprintf('Done!\n');

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
