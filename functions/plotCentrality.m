function p = plotCentrality(G,type)
    %%% Input: Graph G, type of centrality 'type'
    %%% Output: visualize Centrality as size of nodes
    
    rank = centrality(G,type);
    edges = linspace(min(rank),max(rank),7);
    bins = discretize(rank,edges);
    p = plot(G,'Layout','force','EdgeAlpha',0.005,'NodeColor','r')
    p.MarkerSize = bins;
           
end