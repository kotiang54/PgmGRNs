
function A = allcombs(experimentOutcomes)
%
% A function that gives the possible combinations of N input variables
%
% Example: N = 3; # of input variables
%          allcombs(repmat({[1:k]},1, N)); k-logical states
%

    experimentOutcomes = flipud(experimentOutcomes(:));
    n = numel(experimentOutcomes);
    c = cell(1, n);
    [c{:}] = ndgrid(experimentOutcomes{:});
    c = fliplr(c);
    A = cell2mat(cellfun(@(v)v(:), c, 'UniformOutput',false));
end