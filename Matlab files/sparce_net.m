
function new_dag = sparce_net(dag, n)
    % This function creates a sparce network out of the original network
    % It deletes existing edges to create a sparce network 
    % For example we are deleting between 10 and 15 edges.
    

    idx_nz = find(dag); % find the indices of existing edges
    idx_select = datasample(idx_nz, n, 'Replace', false);
    dag(idx_select) = 0;
    new_dag = dag;
    
