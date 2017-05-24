function [ gvecs , qpoints ] = read_gvecs(fname)
% Reads G-vectors corresponding to the chi/epsilon matrix from epsilon.log.
% The result is ordered wrt the epsilon G-space (|q+G|^2) so that the matrix 
% element indexed by e.g. (1,2) naturally corresponds to G=gvecs(1,:), 
% G'=gvecs(2,:) for a given q.
%
% Also returns the corresponding q-points.
%
% FIXME: epsilon.log cannot be from a calculation with skip_epsilon.

fprintf('Reading G-vectors...');

infile = fopen(fname,'r');

qpoints = [];
iq = 0;

while ~feof(infile)
    linestr = fgets(infile);
    parts = textscan(linestr,'%s');
    
    % If the line is not empty and begins with 'q=', it's a q-point
    if (numel(parts{1}) > 1) && (strcmp(parts{1}{1},'q='))
        iq = iq+1;
        qpoints = [qpoints; str2double(parts{1}{2}) str2double(parts{1}{3}) str2double(parts{1}{4})];
    end
    
    % If the line has at least 4 parts and the 4th is 'epsilon', it's the
    % beginning of a set of matrix elements and corresponding G-vectors
    if (numel(parts{1}) >= 4) && (strcmp(parts{1}{4},'epsilon'))
        linestr = fgets(infile); % advance to first line of data
        parts = textscan(linestr,'%f');
        % Read all G-vectors
        gvecs{iq} = [];
        while (numel(parts{1}) == 7) && (~feof(infile))
            gvecs{iq} = [gvecs{iq}; parts{1}(1) parts{1}(2) parts{1}(3)];
            linestr = fgets(infile);
            parts = textscan(linestr,'%f');
        end
        gvecs{iq} = unique(gvecs{iq},'rows','stable'); % 'stable' to maintain ordering
    end
end

fclose(infile);
fprintf('Done.\n');

end