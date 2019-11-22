

function matB = shuffle_data(matA,idx)
% Function to implement shuffling of gene expression data

n = numel(idx);
matB = matA;

if n == 2
   matB(:,[idx(1) idx(2)]) = matA(:,[idx(2) idx(1)]);  
else
    index = randperm(n);
    for i = 1:n
       matB(:,idx(index(i))) = matA(:,idx(i)); 
    end   
end

