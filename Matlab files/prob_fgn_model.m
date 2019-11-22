
function [marginals, count] = prob_fgn_model(data, net_mat,k)

N = size(net_mat,1);
cpt_combs = truth_table(net_mat,k);
regulation_fun = prob_cpd(net_mat, data, k);  % regulation function of the node (CPT-distribution)
full_daG = net_mat + eye(size(net_mat));

% initialize msg updates
var_msg = cell(N,N);
idx_init = find(full_daG == 1);
for i = 1:length(idx_init)
    var_msg{idx_init(i)} = ones(1,k);
%     var_msg{idx_init(i)} = rand(1,k);
end

count = 0;
error = ones(N,1);
tmp1 = zeros(N,1);

while  count < 100 && 1e-4 < max(error)
    % Message Passing updates
    % fact_msg: msgs from factor nodes to variable edges f12 = msg f2 --> x1
    % var_msg: msgs from variable nodes to factor edges  x12 = msg x1 --> f2
    % marginals: posterior distribution of the nodes
    
    fact_msg = f_node_update(regulation_fun,cpt_combs,full_daG,var_msg,k);
    [var_msg, marginals] = v_node_update(full_daG,fact_msg,k);
    
    tmp2 = cell2mat(marginals);
    error = tmp2(:,1) - tmp1;
    tmp1 = tmp2(:,1);
    count = count + 1;
end


