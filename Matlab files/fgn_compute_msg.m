
function msg = fgn_compute_msg(node_function,in_msg)

if isempty(node_function)
    disp('Error at factor node: empty regulation function');
end

msg_snd = in_msg{1};
len = size(in_msg,1);
k_val = size(node_function,2);

for i = 2:len
    msg_snd = kron(msg_snd, in_msg{i});
end
msg_snd = reshape(msg_snd,k_val,[])';
msg_snd = sum(sum(node_function.*msg_snd));

msg = msg_snd;