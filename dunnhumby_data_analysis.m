close all
clear all
clc


addpath functions/
addpath tensor_toolbox/

%load transaction data
filename = 'datasets/transaction_data.csv';
T = readtable(filename);
zeroEntry = find((T.QUANTITY == 0) | (T.SALES_VALUE == 0));
T(zeroEntry,:) = [];

%load product data
product = readtable('datasets/product.csv');
product.PRODUCT_ID=categorical(product.PRODUCT_ID);


T.household_key=categorical(T.household_key);
T.BASKET_ID=categorical(T.BASKET_ID);
T.PRODUCT_ID=categorical(T.PRODUCT_ID);
T.STORE_ID=categorical(T.STORE_ID);

T.B = findgroups(T.BASKET_ID);
T.P = findgroups(T.PRODUCT_ID);
T.H = findgroups(T.household_key);
T.S = findgroups(T.STORE_ID);

aggB = splitapply(@sum,T.QUANTITY,T.B);
plot(sort(aggB, 'descend'));
xlabel('Basket') 
ylabel('Quantity') 

aggP = splitapply(@sum,T.QUANTITY,T.P);
plot(sort(aggP, 'descend'));
xlabel('Product') 
ylabel('Quantity') 

aggH = splitapply(@sum,T.QUANTITY,T.H);
plot(sort(aggH, 'descend'));
xlabel('Household') 
ylabel('Quantity') 

aggS = splitapply(@sum,T.QUANTITY,T.S);
plot(sort(aggS, 'descend'));
xlabel('Product') 
ylabel('Store') 

%get a store in the middle
[out,idx] = sort(aggS, 'descend');

nstore = size(aggS)
store_id = T.STORE_ID(T.S==idx(100))


%get a subset of data of one store and build network
trans = T(T.STORE_ID=='364', :);
[nrow ncol]=size(trans)
trans = join(trans, product);
size(unique(trans.COMMODITY_DESC))

trans.basketID = findgroups(trans.BASKET_ID);
trans.productID = findgroups(trans.PRODUCT_ID);
trans.cmdtyID = findgroups(trans.COMMODITY_DESC);
[n_basket, rows] = size(unique(trans.basketID));
[n_product, rows] = size(unique(trans.productID));
[n_commodity, rows] = size(unique(trans.COMMODITY_DESC));

%build adjacency matrix
A = zeros(n_commodity,n_commodity);

for i = 1 : 10
    basket = trans(trans.basketID == i,:);
    [nrow, ncol] = size(basket);
    for j = 1:(nrow-1)
        for k = (j+1):nrow
            A(basket.cmdtyID(j),basket.cmdtyID(k)) = 1;
            A(basket.cmdtyID(k),basket.cmdtyID(j)) = 1;
        end
    end
end

%save('store364.mat','A')
G = graph(A);
G = max_connected_subgraph(G);
A = adjacency(G);
m = numedges(G);
n = numnodes(G);
deg = degree(G);
edge_density = 2*m/(n*(n-1));
figure, plot(G,'NodeLabel',{});
p = plot(G,'Layout','force','EdgeAlpha',0.005,'NodeColor','r');
names = 1:255;
node_names=char(unique(trans.COMMODITY_DESC));
labelnode(p,names.', names.')

%page rank centrality
pg_ranks = centrality(G,'pagerank');

%degree centrality
deg_ranks = centrality(G,'degree','Importance',G.Edges.Weight);
edges = linspace(min(pg_ranks),max(pg_ranks),7);
bins = discretize(pg_ranks,edges);
p = plot(G,'Layout','force','EdgeAlpha',0.005,'NodeColor','r');
p.MarkerSize = bins;

%build tensor and find higher-order centrality
Tensor = build_triangles_tensor(A,'type','watts_strogatz');
[xarray, resarray] = spectral_cc(A,Tensor,'alpha',0);
c_T = xarray(:,end);
edges = linspace(min(c_T),max(c_T),7);
bins = discretize(c_T,edges);
p = plot(G,'Layout','force','EdgeAlpha',0.005,'NodeColor','r');
p.MarkerSize = bins;

