
function [GRN, genes, exp_data] = gnetwork(grn_type)
% matlab function to create the gene network matrix.
% SOS response network if net=1, and AR network if net=2
% N is the size of the network
% Output: GRN - interaction matrix
%         genes - list of genes
%         exp_data - expression data used
% GRN matrix can be created with any other method known 

switch grn_type
    case 1
        genes = {'lexA','dinI','umuDC','recA','ssb','recF','rpoS','rpoH','rpoD'};  % List of genes
        N = length(genes);
        filename = 'SOS_data.csv';
        GRN = zeros(N,N);
        G1 = 1; G2 = 2; G3 = 3; G4 = 4; G5 = 5; G6 = 6; G7 = 7; G8 = 8; G9 = 9;
        
        GRN(G1, [G2 G5 G6 G9]) = 1;
        GRN(G2, [G4 G6]) = 1;
        GRN(G3, [G1 G2 G5 G6 G9]) = 1;
        GRN(G4, [G1 G3 G5 G6]) = 1;
        GRN(G5, [G2 G6 G9]) = 1;
        GRN(G6, [G7 G9]) = 1;
        GRN(G7, G9) = 1;
        GRN(G8, G9) = 1;
        GRN(G9, [G2 G4]) = 1;
        
    case 2
        genes = {'evgA','ydeO','gadE','gadX','gadA','gadBC','gadW','hdeA',...
            'phoP','crp','rpoS','rpoD','hns','ydeP'};  % List of genes
        N = length(genes);
        filename = 'AR_data.csv';
        GRN = zeros(N,N);
        G1=1; G2=2; G3=3; G4=4; G5=5; G6=6; G7=7; G8=8; G9=9; G10=10;
        G11=11; G12=12; G13=13; G14=14;
        
        GRN(G1, [G2 G3 G14]) = 1;
        GRN(G2, G3) = 1;
        GRN(G3, [G5 G6 G8]) = 1;
        GRN(G4, [G3 G5 G6 G7]) = 1;
        GRN(G7, [G3 G8 G11]) = 1;
        GRN(G9, [G7 G8]) = 1;
        GRN(G10, G11) = 1;
        GRN(G11, G4) = 1;
        GRN(G12, [G5 G6]) = 1;
        GRN(G13, [G1 G11]) = 1;
end

% Expression data
exp_data = csvread(filename, 0,0);

% Network structure and Factor Graph
view(biograph(GRN, genes, 'LayoutType','hierarchical')); % radial, equilibrium, hierarchical
  % or      
G = digraph(GRN,genes);
plot(G,'NodeLabel',genes);
