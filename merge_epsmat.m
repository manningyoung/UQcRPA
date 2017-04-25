function [epsmat_full , nmtx_full] = merge_epsmat(epsmat,eps0mat,nmtx,nmtx0)
% Constructs the full inverse dielectric matrix from eps0mat and epsmat.
% Also removes the padding zeros from epsmat and returns nmtx for the 
% full matrix.

% The first cell is epsilon^{-1}(q=q0)
epsmat_full{1} = reshape(eps0mat,[nmtx0 nmtx0]);

% And the rest are epsilon^{-1}(q)
for iq = 1:length(nmtx)
    epsmat_full{iq+1} = reshape(epsmat(1,1:nmtx(iq),1:nmtx(iq),1,1,iq),[nmtx(iq) nmtx(iq)]);
end

nmtx_full = [nmtx0; nmtx];

end