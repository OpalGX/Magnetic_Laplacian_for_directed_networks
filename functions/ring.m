function [G, A] = ring(k, n, p_in, p_out)
%generatl a k ring with n*k nodes as benchmark for directed graph

A = zeros(k*n,k*n);
for i = 1:n*k
    r = ceil(i/n);
    for j = 1:n*k
        c = ceil(j/n);
        if i == j %diagonal entries are zeros
            A(i,j) = 0;
        elseif c == r %same cluster
            A(i,j) = (rand(1) <= p_in);
        elseif c == r+1
            A(i,j) = (rand(1) <= p_out);  
        elseif c == mod(r+1,k)
            A(i,j) = (rand(1) <= p_out);
        %elseif r == mod(c+1,k)
        %    A(i,j) = (rand(1) <= 0.05);
        else
            A(i,j) = 0;
        end
    end
end
%idx = randperm(k*n);
%A = A(idx,idx);

G = digraph(A);
G.Nodes.L = reshape(repmat([1:k], n, 1), [], 1); %categories of nodes
%G = reordernodes(G,randperm(k*n)); %randomize matrix
end
