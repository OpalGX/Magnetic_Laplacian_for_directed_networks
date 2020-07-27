function [G, eig, p, V, D] = meigenmaps(G,g)

    W = adjacency(G);
    % Symmetric weights.
    Ws = (W+W.')/2;

    %Edge flow.
    A = W-W.';
    Gs = graph(Ws);

    deg = sum(Ws,2); Deg = diag(deg); % Degree matrix.?
    Tg = exp(1) .^(2*pi*1i*g*A.');  % Transporter
    Lg = Deg - Ws.*Tg; % Magnetic Laplacian.
    Lng = inv(sqrtm(Deg))*Lg*inv(sqrtm(Deg));% Normalized Laplacian.
    [V,D] = eigs(Lng,size(W,1),'smallestabs'); % Eigenvectors ranging from smallest
    Phi = mod(angle(V),2*pi); %Phases.
 
    G.Nodes.phase0 = Phi(:,1);
    G.Nodes.phase1 = Phi(:,2);
    
    %p = plot(G,'XData', G.Nodes.phase0 ,'YData',G.Nodes.phase1 , 'NodeLabel', G.Nodes.Name);
    p = plot(G,'XData', cos(G.Nodes.phase0) ,'YData', sin(G.Nodes.phase0));
    p.Marker = 's';
    p.NodeColor = 'r';
    xlabel('phase(\phi_0)');
    ylabel('phase(\phi_1)');
    
    eig = round(D(1,1),4);
    

end
