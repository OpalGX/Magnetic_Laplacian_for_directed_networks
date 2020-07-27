close all
clear all
clc

addpath functions/
addpath tensor_toolbox/

tic

%load transaction data
filename = 'datasets/transaction_data.csv';
T = readtable(filename);
zeroEntry = find((T.QUANTITY == 0) | (T.SALES_VALUE == 0)); %rem
T(zeroEntry,:) = [];

%load product data
product = readtable('datasets/product.csv');

%aggregate quantity by store
T.S = findgroups(T.STORE_ID);
aggS = splitapply(@sum,T.QUANTITY,T.S);




%rank store by quantity
[out,idx] = sort(aggS, 'descend');
nstore = size(aggS);
store_id = T.STORE_ID(T.S==idx(100));


%get a subset of data of one store and build network
store_data = T(T.STORE_ID==364, :);
store_data = join(store_data, product);
all_trans = join(T, product);



%extract high frequency commodities from all stores
%aggregate quantity by commodity
all_trans.C = findgroups(all_trans.COMMODITY_DESC);
aggC = splitapply(@sum,all_trans.QUANTITY,all_trans.C);

%aggregate quantity by SUB_COMMODITY_DESC
all_trans.Sub = findgroups(all_trans.SUB_COMMODITY_DESC);
aggSub = splitapply(@sum,all_trans.QUANTITY,all_trans.Sub);


%filter basket with the high frequency commodities
%rank store by quantity
[out,idx] = sort(aggC, 'descend');
nCmdty = size(aggC);
topCmdty = idx(1:ceil(10*nCmdty/100));
top_cmdty_trans = all_trans(ismember(all_trans.C, topCmdty), :);

save('processed_data/trans_data','all_trans', 'top_cmdty_trans','store_data')

toc

