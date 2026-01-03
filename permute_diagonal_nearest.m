function [D2_permuted, perm_idx] = permute_diagonal_nearest(D1, D2)
    % PERMUTE_DIAGONAL_NEAREST Permute diagonal matrix D2 based on nearest values in D1
    %
    % Inputs:
    %   D1 - First diagonal matrix (n x n)
    %   D2 - Second diagonal matrix (n x n) to be permuted
    %
    % Outputs:
    %   D2_permuted - Permuted version of D2
    %   perm_idx - Permutation indices used
    
    % Extract diagonal elements
    d1 = diag(D1);
    d2 = diag(D2);
    n = length(d1);
    
    % Initialize permutation index
    perm_idx = zeros(n, 1);
    used = false(n, 1);
    
    % For each element in d1, find the nearest unused element in d2
    for i = 1:n
        % Compute distances to all unused elements in d2
        distances = abs(d2 - d1(i));
        distances(used) = inf;  % Mark used indices as unavailable
        
        % Find the nearest element
        [~, nearest_idx] = min(distances);
        
        % Record the permutation
        perm_idx(i) = nearest_idx;
        used(nearest_idx) = true;
    end
    
    % Apply permutation to D2
    d2_permuted = d2(perm_idx);
    D2_permuted = diag(d2_permuted);
end