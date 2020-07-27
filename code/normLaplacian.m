function [G, p] = normLaplacian(G)


    W = adjacency(G);
    Ws = (W+W.')/2;
    deg = outdegree(G); D = diag(deg); % Degree matrix.?
    L = D - Ws; % Magnetic Laplacian.
    Ln = inv(sqrtm(D))*L*inv(sqrtm(D));% Normalized Laplacian.
    [V,D] = eig(Ln); % Eigenvectors.

    G.Nodes.Eig0 = V(:,1);
    G.Nodes.Eig1 = V(:,2);
    G.Nodes.Eig2 = V(:,3);

    %p = plot(G,'XData', G.Nodes.Eig0 ,'YData',G.Nodes.Eig1 , 'ZData',G.Nodes.Eig2, 'NodeLabel', G.Nodes.Name);
    p = plot(G,'XData', G.Nodes.Eig0 ,'YData',G.Nodes.Eig1 ,  'NodeLabel', G.Nodes.Name);

    p.Marker = 's';
    p.NodeColor = 'r';
    xlabel('Eig_0');
    ylabel('Eig_1');
    zlabel('Eig_2');


end
