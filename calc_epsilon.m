function [ epsmat ] = calc_epsilon(chi0mat,chimat,chi0mat_a,chimat_a,vcoul0,vcoul,nmtx0,nmtx)
% Constructs the contstrained inverse dielectric matrix from the full
% polarisability P and that of the active subspace P_a according to
%   epsilon_r^{-1} = 1 - v*(P-P_a)
% P and P_a must be the same size (nmtx).

nmtx = [nmtx0; nmtx]; % append nmtx(q0) to the list

for iq = 1:length(nmtx)
    for ig = 1:nmtx(iq)
        for igp = 1:nmtx(iq)
            if iq==1 % handle q0 separately
                epsmat{iq}(ig,igp) = double(ig==igp) - vcoul0(ig,iq)*(chi0mat(1,ig,igp,1,1,iq)-chi0mat_a(1,ig,igp,1,1,iq));
            else
                epsmat{iq}(ig,igp) = double(ig==igp) - vcoul(ig,iq)*(chimat(1,ig,igp,1,1,iq)-chimat_a(1,ig,igp,1,1,iq));
            end
        end
    end
    epsmat{iq} = inv(epsmat{iq});
end
end