function [ unk ] = read_bloch(fname)
% Reads the periodic part of the Bloch states u_{n,k} from UNKp.s, where
% p is the k-point index and s is the spin (= 1 or 2).

fprintf('Reading %s...',fname);

infile = fopen(fname,'r');

% Parse header
linestr = fgets(infile);
parts = textscan(linestr, '%f');
ngx = parts{1}(1); % number of grid points in each direction
ngy = parts{1}(2);
ngz = parts{1}(3);
ik = parts{1}(4); % k-point index
num_band = parts{1}(5); % total number of bands

iband = 1; 
ir = 1;

while ~feof(infile)
    linestr = fgets(infile);
    parts = textscan(linestr, '%f');
    if numel(parts{1}) == 2 % is a value
        unk{iband}(ir) = complex(parts{1}(1),parts{1}(2));
        ir = ir+1;
        if ir > ngx*ngy*ngz % goto new band
            iband = iband+1;
            ir = 1;
        end
        if iband > num_band+1
            error('Number of bands exceeded while reading %s. Possible formatting error.',fname)
        end
    end
end

fprintf('Done.\n')

end