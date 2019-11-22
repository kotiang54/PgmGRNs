
function msg = normalization(msg_in)

    % Normalization of the messages
    
    sum_msg = sum(msg_in);
    if sum_msg == 0
        msg = ones(size(msg_in)); % division of zero error
    else
        msg = msg_in/sum_msg;
    end
      
end