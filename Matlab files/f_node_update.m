
function fact_msg = f_node_update(regulation_fun,cpt_combs,daG_mat,input_msg,k)

% inputs: regulation_fun: gene regulation function
%         cpt_combs: the combinatoria truth table
%         idx: node ids
%         daG_mat: full adjacency matrix of the network
%         input_msg: input messages use in computation
%         k: k-dim message space, logical levels

% output: nxn updated factor node messages (fact_msg)

N = size(daG_mat,1);
fact_msg = cell(N,N);

% Factor node message passing updates
for i = 1:N
    % Edges from factor nodes
    idx = find(daG_mat(:,i) == 1);
    msg_snd = zeros(1,k);
    if length(idx) == 1   % leaf factor node
        msg_snd = regulation_fun{i};
        fact_msg{idx,i} = msg_snd;      % store the msg somewhere!!!
    else
        var = find(idx == i);
        new_idx = idx;
        new_idx(var) = [];
        for j = 1:length(idx)
            % j is the destination of message
            % collect all messages except from j
            msg_in = {};
            for m = 1:length(idx)
                if j == m
                    continue;
                end
                msg_in = [msg_in; input_msg{idx(m), i}];
            end
            % calculate the message HERE !!!!!
            node_potential = [];
            if idx(j) == i
                for n = 1:k
                    node_potential = reshape(regulation_fun{i}(:,n),k,[])';
                    msg_snd(n) = fgn_compute_msg(node_potential, msg_in);
                end
            else
                tmp = j;    % temporary store of index j
                if j > var
                    tmp = j-1;
                end
                var_idx = find(new_idx == new_idx(tmp));  % identity of node sending msg to
                for p = 1:k
                    for q = 1:size(cpt_combs{i},1)
                        if cpt_combs{i}(q,var_idx) == p
                            node_potential = [node_potential; regulation_fun{i}(q,:)];
                        end
                    end
                    msg_snd(p) = fgn_compute_msg(node_potential, msg_in);
                    node_potential = [];
                end
            end
            fact_msg{idx(j),i} = normalization(msg_snd);  % Normalization of the messages
            msg_snd = [];
        end
    end
end








