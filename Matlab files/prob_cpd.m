
%
% Conditional Probability Distribution calculation
% inputs: network, net; and the number of logical states, k. 
% data, represents the discretized gene expression data in k logical states
% The cpd table looks like: eg for k = 2, with 2 inputs (G1|G2;G3)
%      G2 G3     G1     
%           |  1    2   |
%     ------|-----------|---
%      1  1 | 0.7  0.3  |
%      1  2 | 0.25 0.75 |
%      2  1 | 0.01 0.99 |
%      2  2 | 1.0  0.0  |
%     ------|-----------|---     
%

function cpd = prob_cpd(net, data, k)

if size(net,1) ~= size(net,2)
    error('Invalid network');
else
    n = size(net,1);
    cpd = cell(1, n);
    
    for i = 1:n
        idx = find(net(:,i) == 1);
        r = numel(idx);
        
        % Matrix of all possible combinations given n inputs and k
        % logical states
        pos_combs = allcombs(repmat({1:k}, 1, r+1));
        data_select = [data(:,idx) data(:,i)];
        y = k^r;    % total number of rows in the cpd tables
        joint_pr_dist = zeros(y, k);
        index = 0;
        
        for t = 1:y           
            new_combs = pos_combs(1+index :k*t,:);     % select the combinations to apply every k logical states
            if r == 0   % r=0 represents node(s) without external influence
                count_evid = size(data_select,1);  % Count the total number of occurances
                for j = 1:k
                    count_indv = numel(find(data_select == j));
                    joint_pr_dist(t,j) = count_indv./count_evid;
                end            
            elseif r == 1  % r=1, node(s) with only 1 external influence
                joint_pr_dist(t,:) = calculate_prob(data_select, new_combs, k);               
            else
                joint_pr_dist(t,:) = calculate_prob(data_select, new_combs, k);
            end
            index = index + k;
        end
        cpd{i} = joint_pr_dist./sum(joint_pr_dist,2);
        cpd{i}(isnan(cpd{i})) = 0;
    end
end
end

