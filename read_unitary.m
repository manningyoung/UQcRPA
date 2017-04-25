function [ unitary_matrix , kpoints ] = read_unitary(fname)
% seedname_u.mat contains the matrix elements of U^(k) in column-major
% order for each k-point.
%
% Requires write_u_matrices=.true. in win.

infile = fopen(fname,'r');

linestr = fgets(infile); % date info

% read matrix info
linestr = fgets(infile);
parts = textscan(linestr, '%f');
num_kpts = parts{1}(1); % number of k-points
num_wann = parts{1}(2); % number of wannier functions

kpoints = []; % might as well read the k-points while we're at it

linestr = fgets(infile); % advance
ik = 0;
while ischar(linestr) % until EOF
    parts = textscan(linestr, '%f');
    if numel(parts{1}) == 3 % is a k-point
        kpoints = [kpoints; parts{:}'];
        row = 1; col = 1; ik = ik+1; % reset indices
    elseif numel(parts{1}) == 2 % is a matrix element
        unitary_matrix{ik}(row,col) = complex(parts{1}(1),parts{1}(2));
        row = row+1;
        if row > num_wann % goto next col
            col = col+1; row = 1;
        end
    end
    linestr = fgets(infile); % advance
end

end