
% Main file to implement network perturbation method (shuffled nodes)
% Handles both Ecoli SOS response network , and Acid Resistance (AR)
% regulatory network
% Project: " A PGM for system-wide analysis of GRNs"
% by S. Kotiang and A. Eslami

clear all;
clc;
warning off

prompt = 'Type 1 for SOS and 2 for AR -> ';
net =  sscanf(input(prompt, 's'), '%d');
fprintf('\n');

prompt = 'Select the number of clusters (2 or 3) k -> ';
k =  sscanf(input(prompt, 's'), '%d');
fprintf('\n');

% Network structure and Factor Graph
[Ecoli_dag, genes, exp_data] = gnetwork(net); % create the GRN interaction matrix
full_daG = Ecoli_dag + eye(size(Ecoli_dag)); % full network matrix
[data_len, N] = size(exp_data);

[dEcoli, class_proportions] = discretization(exp_data,k, net); % discretization of the expression data
X = class_proportions(:);

% update the network
fgn_edges = nnz(Ecoli_dag);
n_rand = 0:5;
len = length(n_rand);
corr_r = zeros(1,len);
iterations = zeros(1,len);
idx1 = find(Ecoli_dag);  % position index of ones in the matrix
idx0 = find(full_daG == 0); % position index of zeros in the matrix
diag_idx = logical(eye(size(Ecoli_dag)));
index_diag = find(diag_idx);

for i = 1:len
    rand_attempts = 300; % or: fgn_edges * 10
    if n_rand(i) == 0
        [marginals, count] = prob_fgn_model(dEcoli,Ecoli_dag,k);
        beliefs = cell2mat(marginals);
        Y = beliefs(:);
        r = corrcoef(X,Y,'alpha',0.01);
        corr_r(i) = r(1,2);
        iterations(i) = count;
    else
        attemps = nchoosek(fgn_edges,n_rand(i));
        if attemps < rand_attempts
            rand_attempts = attemps;
        end
        iter = zeros(1,rand_attempts);
        coeffs = zeros(1,rand_attempts);
        
        for j = 1:rand_attempts
            random_net = Ecoli_dag;
            pick0 = randsample(idx1, n_rand(i)); % randomly select edges to swap
            pick1 = randsample(idx0, n_rand(i)); % randomly select edges to swap
            random_net(pick0) = 0;
            random_net(pick1) = 1;
            
            % Network computation
            [marginals, count] = prob_fgn_model(dEcoli,random_net,k);
            beliefs = cell2mat(marginals);
            Y = beliefs(:);
            r = corrcoef(X,Y,'alpha',0.01);
            coeffs(j) = r(1,2);
            iter(j) = count;
        end
        corr_r(i) = mean(coeffs);
        iterations(i) = mean(iter);
    end
end

figure(2);
plot(0:5,corr_r, 'o-b');

% Store the correlation rho values 
% corr_r_2 = [0.9883 0.9742 0.9304 0.8978 0.8880 0.8969];  % SOS network
% corr_r_3 = [0.9611 0.9555 0.9394 0.9258 0.9354 0.9153];  % SOS network
% 
% corr_r_2 = [0.9298 0.9100 0.8831 0.8710 0.8517 0.8493];  % AR network
% corr_r_3 = [0.9398 0.9300 0.8869 0.8763 0.8618 0.8573];  % AR network
% % 
% plot(0:5,corr_r_2, 'o-b', 0:5, corr_r_3, 'o-r')
% xlabel('n-randomized edges');  ylabel('Average correlation coeff.');
% legend('2-states','3-states')
 


