
function [dct_data, class_proportions, gmfit]  = discretization(data,k, grn_type)

% function takes in continuous expression data and output discrete data
% Output: dct_data - discretized data
%         class_proportions - distribution of the discretized data per gene
%         node
%         gmfit - learned GMM function

options = statset('MaxIter',1000);
gmfit = {};
rndn1 = [1236 981 272 14 1];    % random numbers for reproducibility
rndn2 = [276 196 10 9 2 18 3 24];
[data_len, N] = size(data);
dct_data = zeros(data_len, N);
class_proportions = zeros(1,k); % discretized class proportions

for i = 1: N
    switch grn_type
        case 1
            if k >= 2 && k <= 6
                rng(rndn1(k-1));
            end
        case 2
            if k >= 2 && k <= 9
                rng(rndn2(k-1));
            end           
        otherwise rng default;
    end
%     gmfit{i} = fitgmdist(data(:,i),k, 'start', 'plus','Options', options, 'Replicates',10,'RegularizationValue',0.01);
    gmfit{i} = fitgmdist(data(:,i),k, 'start', 'plus','Options', options, 'Replicates',10);
    dct_data(:,i) = cluster(gmfit{i}, data(:,i));
    for j = 1:k
        class_proportions(i,j) = numel(find(dct_data(:,i)==j))/data_len;
    end
end