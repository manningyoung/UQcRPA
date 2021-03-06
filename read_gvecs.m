function [ gvecs , qpoints ] = read_gvecs(fname)
% Reads G-vectors corresponding to the chi/epsilon matrix from epsilon.log.
% The result is ordered wrt the epsilon G-space (|q+G|^2) so that the matrix 
% element indexed by e.g. (1,2) naturally corresponds to G=gvecs(1,:), 
% G'=gvecs(2,:) for a given q.
%
% Also returns the corresponding q-points.

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
    
    % If the line has at least 5 parts and the 5th is 'chi(g,gp)', it's the
    % beginning of a set of matrix elements and corresponding G-vectors
    if (numel(parts{1}) >= 5) && (strcmp(parts{1}{5},'chi(g,gp)'))
        gvecs{iq} = [];
        while (numel(parts{1}) > 0) && (~feof(infile))
            linestr = fgets(infile); % advance to first line of data
            parts = textscan(linestr,'%f');
            if numel(parts{1} == 9) % check we're still reading G-vectors
                gvecs{iq} = [gvecs{iq}; parts{1}(1) parts{1}(2) parts{1}(3)];
            end
        end
        gvecs{iq} = unique(gvecs{iq},'rows','stable'); % 'stable' to maintain ordering
    end
    
end

fclose(infile);
fprintf('Done.\n');

end