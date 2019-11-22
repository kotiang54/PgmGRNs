
function [msg, belief] = v_node_update(daG_mat,input_msg,k)

% inputs:
%         daG_mat: full adjacency matrix of the network
%         input_msg: input messages use in computation
%         k: k-dim message space, logical levels

% output: nxn updated factor node messages (var_msg)
%         nx1 posterior marginals (belief) of the node(var_msg)
% variable node simply multiply the product of all incoming
% messages except the i-th node being sent to

N = size(daG_mat,1);
msg = cell(N,N);
belief = cell(N,1);

% Factor node message passing updates
for i = 1:N
    % Edges from variable nodes
    idx = find(daG_mat(i,:) == 1);   
    if length(idx) == 1   % leaf variable node
        msg_snd = ones(1,k);
        msg{i,idx} = normalization(msg_snd);  % store the msg somewhere!!!
    else
        for j = 1:length(idx)
            % j is the destination of message
            % collect all messages except from j
            msg_in = {};
            msg_snd = ones(1,k);
            for m = 1:length(idx)
                if j == m
                    continue;
                end
                msg_in = [msg_in; input_msg{i,idx(m)}];
            end
            len = size(msg_in,1);
            for m = 1:len
                msg_snd = msg_snd.*msg_in{m};
            end
            msg{i,idx(j)} = normalization(msg_snd);
        end      
    end 
    
    % calculate the node belief: collect all messages
    % variable node simply multiply the product of all incoming
    
    % how to treat isolated nodes.
    if nnz(daG_mat(i,:))== 1 && nnz(daG_mat(:,i))== 1 % case of isolated nodes
        belief{i,1} = zeros(1,k);
    else
        msg_bel = ones(1,k);
        for t = 1:length(idx)
            msg_bel = msg_bel.*input_msg{i,idx(t)};
        end
        belief{i,1} = normalization(msg_bel);
    end
  
end







