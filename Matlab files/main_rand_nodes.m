
% Main file to implement network perturbation method (shuffled data)
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

% initialize FGN msg updates
var_msg = cell(N,N);  % variable node messages
idx_init = find(full_daG == 1);
for i = 1:length(idx_init)
    var_msg{idx_init(i)} = ones(1,k);
    %     var_msg{idx_init(i)} = normalization(rand(1,k));
end

% update the network
n_rand = 0:N;
len = length(n_rand);
corr_r = zeros(1,len-1);
iterations = zeros(1,len-1);

for i = 1:len-2
    rand_attempts = 100;
    if n_rand(i) == 0  % true network without perturbation
        [marginals, count] = prob_fgn_model(dEcoli,Ecoli_dag,k);
        beliefs = cell2mat(marginals);
        Y = beliefs(:);
        r = corrcoef(X,Y,'alpha',0.01);
        corr_r(i) = r(1,2);
        iterations(i) = count;
    else  % perturbation of the network nodes
        attempts = nchoosek(N,n_rand(i));
        if attempts < rand_attempts
            rand_attempts = attempts;
        end
        iter = zeros(1,rand_attempts);
        coeffs = zeros(1,rand_attempts);
        
        for j = 1:rand_attempts
            idx = randperm(N,n_rand(i)+1);
            data = shuffle_data(dEcoli,idx);
            
            % Network computation
            [marginals, count] = prob_fgn_model(data,Ecoli_dag,k);
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

% Plots
figure(2);
axis square;
axis([0 1 0 1]);
plot(0:N-1,corr_r, 'o-b');

% STORE CORRELATION VALUES OF 2- AND 3_STATES AND PLOT ON SAME GRAPH

% figure(3)
% corr_r_2 = [0.9883 0.9111 0.8189 0.6726 0.6822 0.4909 0.4125 0.3827 0];
% corr_r_3 = [0.9611 0.6087 0.7747 0.7331 0.6344 0.5790 0.5288 0.4669 0];
% plot(0:N-1,corr_r_2, 'o-b', 0:N-1, corr_r_3, 'o-r')
% legend('2-states','3-states')

