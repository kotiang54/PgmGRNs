
function possible_comb = truth_table(network,k)
% Function to create a truth table for CPT computation

N = size(network,1);
possible_comb = cell(1,N);
for i = 1:N
    idx = find(network(:,i)==1);
    r = numel(idx);
    if r == 0
        possible_comb{i} = [];
    else
        possible_comb{i} = allcombs(repmat({1:k}, 1, r));
    end
end