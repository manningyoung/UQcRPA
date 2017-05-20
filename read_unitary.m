function [ unitary_matrix , kpoints ] = read_unitary(fname)
% Reads the unitary transformation matrix U^(k) from seedname_u.mat.
% Matrix elements are read in column-major order for each k-point.
%
% Also returns the corresponding k-points.

fprintf('Reading unitary MLWF transformation matrix...');

infile = fopen(fname,'r');

% Parse header
linestr = fgets(infile); % date info
linestr = fgets(infile);
parts = textscan(linestr, '%f');
num_kpts = parts{1}(1); % number of k-points
num_wann = parts{1}(2); % number of wannier functions

kpoints = []; % read the k-points while we're at it
ik = 0;

while ~feof(infile)
    linestr = fgets(infile);
    parts = textscan(linestr, '%f');
    if numel(parts{1}) == 3 % is a k-point
        kpoints = [kpoints; parts{:}'];
        row = 1; 
        col = 1; 
        ik = ik+1;
    elseif numel(parts{1}) == 2 % is a matrix element
        unitary_matrix{ik}(row,col) = complex(parts{1}(1),parts{1}(2));
        row = row+1;
        if row > num_wann % goto next col
            col = col+1; 
            row = 1;
        end
    end
end

fclose(infile);
fprintf('Done.\n');

end