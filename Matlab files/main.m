
% Main file for the gene steady-state distribution simulation
% Handles both Ecoli SOS response network , and Acid Resistance (AR)
% regulatory network.
% Can be adapted to any GRN given gene expression data.
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
pos_comb = truth_table(Ecoli_dag,k);
node_fun = prob_cpd(Ecoli_dag, dEcoli, k);  % regulation function of the node (CPT-distribution)

% initialize FGN msg updates
var_msg = cell(N,N);  % variable node messages
idx_init = find(full_daG == 1);
for i = 1:length(idx_init)
    var_msg{idx_init(i)} = ones(1,k);
    %     var_msg{idx_init(i)} = normalization(rand(1,k));
end

% update the network
count = 0;
error = ones(N,1);
tmp1 = zeros(N,1);
corr_r = [];

tic
while  count < 100 && 1e-4 < max(error)
    % Message Passing updates
    % fact_msg: msgs from factor nodes to variable edges f12 = msg f2 --> x1
    % var_msg: msgs from variable nodes to factor edges  x12 = msg x1 --> f2
    % marginals: posterior distribution of the nodes
    
    fact_msg = f_node_update(node_fun,pos_comb,full_daG,var_msg,k);
    [var_msg, marginals] = v_node_update(full_daG,fact_msg,k);
    
    tmp2 = cell2mat(marginals);
    error = tmp2(:,1) - tmp1;
    tmp1 = tmp2(:,1);
    count = count + 1;
    
    % Keep track of correlation coefficients
    % Uncomment the lines below for plots
    
%     beliefs = cell2mat(marginals);
%     Y = beliefs(:);
%     [r, p] = corrcoef(X,Y,'alpha',0.01);
%     corr_r = [corr_r; r(1,2)];
end
toc

beliefs = cell2mat(marginals);
Y = beliefs(:);
[r, p] = corrcoef(X,Y,'alpha',0.01);
corr_coef = r(1,2);
p_value = p(1,2)
iterations = count
genes = genes';
Tab = table(genes, beliefs, class_proportions);
disp(Tab)

% Correlation plots
figure(2);
axis square;
axis([0 1 0 1]);
scatter(X,Y, 'ob');
lsline;
ylabel('loopy belief');
xlabel('node proportions');

str = ['r= ',num2str(corr_coef)];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
set(T, 'verticalalignment', 'top', 'horizontalalignment', 'left');

% Iteration plots of the correlation coeffs
% How rhos improve with increasing iteration
% figure(3);
% axis square;
% axis([0 1 0 1]);
% plot(1:length(corr_r), corr_r,'r-o')
% ylabel('Correlation coeff. \rho ')
% xlabel('iterations')
