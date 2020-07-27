function [score] = calculate_clustering_precision(A, k, n)

for distance = 1:k-1;
    %clockwise
    for i = 1:k
        submatrix = A((i-1)*n+1:i*n, (mod((i+distance-1)*n+1, n*k):mod((i+distance-1)*n+1, n*k)+n-1));
        n_ones_clockwise(i)= sum(submatrix, 'all');
    end


    total(distance) = sum(n_ones_clockwise);

end

score = max(total)/sum(A,'all');

end
