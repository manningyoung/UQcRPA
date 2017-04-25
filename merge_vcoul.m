function [ vcoul_full ] = merge_vcoul(vcoul,vcoul0,nmtx)
% Constructs the bare coulomb interaction from vcoul(q=q0) and vcoul(q)
% and removes padding zeros.

vcoul_full{1} = vcoul0;

for iq = 1:length(nmtx)
    vcoul_full{iq+1} = vcoul(1:nmtx(iq),iq);
end

end