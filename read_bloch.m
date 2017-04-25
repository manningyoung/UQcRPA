function [ unk ] = read_bloch(fname)
% Reads the wannier90 output file UNKp.s containing the periodic part of the
% Bloch states u_{nk} indexed by k-point p (=1:num_kpts) and spin s (=1 or 2).
%
% Requires write_unk=.true. and wvfn_formatted=.true. in pw2wan.

infile = fopen(fname,'r');

% File header info
linestr = fgets(infile);
parts = textscan(linestr, '%f');
ngx = parts{1}(1); % number of grid points in each direction
ngy = parts{1}(2);
ngz = parts{1}(3);
ik = parts{1}(4); % k-point index
num_band = parts{1}(5); % total number of bands

% Wavefunctions
iband = 1; ir = 1;
while ischar(linestr) % until EOF
    parts = textscan(linestr, '%f');
    if numel(parts{1}) == 2 % is a value
        unk{iband}(ir) = complex(parts{1}(1),parts{1}(2));
        ir = ir+1; % increment spatial index
        if ir > ngx*ngy*ngz % spatial grid exceeded, goto new band
            iband = iband+1;
            ir = 1;
        end
        if iband > num_band+1
            error('Number of bands exceeded while reading %s. Possible file formatting error.',fname)
        end
    end
    linestr = fgets(infile); % advance
end

end