
function new_dag = dense_net(dag, n, diag_ind)
    % This function creates a dense network out of the original network
    % It adds edges to create a dense network 
    % For example we randomly add between 30 and 40 edges.
    % Inputs:
             % dag - matrix of the network
             % diag_ind - indices of the diagonal matrix
    % output: new_dag - dense network created        
       
    
    idx_nz = find(~dag); % find the indices of non-existing edges
    nums = ismember(idx_nz, diag_ind);
    idx_nz(nums) = [];   % remove the indices of diagonal elements
    idx_select = datasample(idx_nz, n, 'Replace', false);
    dag(idx_select) = 1;
    new_dag = dag;