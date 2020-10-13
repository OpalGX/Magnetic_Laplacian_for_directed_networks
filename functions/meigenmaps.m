function [G, eig, p, Phi, D] = meigenmaps(G,g)
%G plot the nodes against the eigenvalue
%eig is the top eigenvector
%p is the network plot
%Phi are the phases corresponding to the top eigenvector
%D is the top eigenvalue

    W = adjacency(G);
    % Symmetric weights.
    Ws = (W+W.')/2;

    %Edge flow.
    A = W-W.';
    Gs = graph(Ws);

    deg = sum(Ws,2); Deg = diag(deg); % Degree matrix.?
    Tg = exp(1) .^(2*pi*1i*g*A.');  % Transporter
    Lg = Deg - Ws.*Tg; % Magnetic Laplacian.
    %Lg = inv(sqrtm(Deg))*Lg*inv(sqrtm(Deg));% Normalized Laplacian.
    [V,D] = eigs(Lg,size(W,1),'smallestabs'); % All eigenvectors ranging from smallest eigenvalue to largest eigenvalue
    D = diag(D);
    Phi = angle(V); %phases corresponding to all eigenvectors
 
    G.Nodes.phase0 = mod(Phi(:,1), 2*pi);%phases corresponding to the top eigenvector
    G.Nodes.phase1 = mod(Phi(:,2), 2*pi);
    
    %p = plot(G,'XData', G.Nodes.phase0 ,'YData',G.Nodes.phase1 , 'NodeLabel', G.Nodes.Name);
    p = plot(G,'XData', cos(G.Nodes.phase0) ,'YData', sin(G.Nodes.phase0)); %plot on polar coordinate
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
    eig = round(D(1,1),4); %top eigenvalue
end
