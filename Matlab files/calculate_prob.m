
%%
%  A function to calculate the probability of an occurance
%  The inputs: data, prior evidence probability, k logical states 
%  and the combinations considered. 
%
%%
function [prob_cal] = calculate_prob(data, combs, k)
    
    evidence = size(data,1);
    prob_cal = zeros(1, k);
    for j = 1:k
        if evidence == 0
            prob_cal(1,j) = 0;
        else
            % Count the total number of occurances
            count_data = sum(ismember(data, combs(j,:),'rows'));
            prob_cal(1,j) = count_data./evidence;
        end
    end

end