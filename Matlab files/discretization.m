
function [dct_data, class_proportions, gmfit]  = discretization(data,k, grn_type)

% function takes in continuous expression data and output discrete data
% Output: dct_data - discretized data
%         class_proportions - distribution of the discretized data per gene
%         node
%         gmfit - learned GMM function

options = statset('MaxIter',1000);
gmfit = {};
[data_len, N] = size(data);
dct_data = zeros(data_len, N);
class_proportions = zeros(1,k); % discretized class proportions

for i = 1: N
    switch grn_type
        case 1
            if k == 2, rng(1236), elseif k == 3, rng(981), end
        case 2
            if k == 2, rng(276), elseif k == 3, rng(196), end
            
        otherwise rng default;
    end
    gmfit{i} = fitgmdist(data(:,i),k, 'start', 'plus','Options', options, 'Replicates',10);
    dct_data(:,i) = cluster(gmfit{i}, data(:,i));
    for j = 1:k
        class_proportions(i,j) = numel(find(dct_data(:,i)==j))/data_len;
    end
end