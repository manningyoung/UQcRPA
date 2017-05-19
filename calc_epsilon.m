function [epsmat, vcoul_full] = calc_epsilon(chi0mat,chimat,chi0mat_a,chimat_a,nmtx0,nmtx,vcoul0,vcoul)
% Constructs the constrained inverse dielectric matrix from the full
% polarisability P and active polarisability P_a according to
%   epsilon^{-1} = [1-v(P-P_a)]^{-1}.
%
% Returns a cell array epsmat containing an nmtx(q)-by-nmtx(q) matrix for each
% q (including q0), with no padding. Also returns corresponding vcoul.

% Handle q0 and store it in the first cell of epsmat
for ig = 1:nmtx0
    vcoul_full{1}(ig) = vcoul0(ig,1);
    for igp = 1:nmtx0
        epsmat{1}(ig,igp) = double(ig==igp) - vcoul0(ig,1)*(chi0mat(1,ig,igp)-chi0mat_a(1,ig,igp));
    end
end

% And the remaining q
for iq = 1:length(nmtx)
    for ig = 1:nmtx(iq)
        vcoul_full{iq+1}(ig) = vcoul(ig,iq);
        for igp = 1:nmtx(iq)
            epsmat{iq+1}(ig,igp) = double(ig==igp) - vcoul(ig,iq)*(chimat(1,ig,igp,1,1,iq)-chimat_a(1,ig,igp,1,1,iq));
        end
    end
end

% Invert
for iq = 1:length(nmtx)+1
    % Check the condition
    if (rcond(epsmat{iq}) < 10*eps) % eps = machine epsilon
        warning('Badly conditioned epsilon matrix. Results may be inaccurate.');
    end
    epsmat{iq} = inv(epsmat{iq});
end
end