function phi = grid2d_to_phi(grid2d)

dim1 = size(grid2d, 1); dim2 = size(grid2d, 2);

% Vectorize the sampling grid
vect_sampling = reshape(grid2d, [dim1*dim2, 1]);

N = sum(vect_sampling(:));       % the number of the sparse sampling points

p_sample = find(vect_sampling == 1 );

% Under-sampling matrix
phi = zeros(N, dim1*dim2);    % Sensing matrix

for i = 1:N
    phi(i, p_sample(i)) = 1;
end

% Visualizing the under-sampling matrix
% figure; imagesc(phi)
% figure; imagesc(grid2d)

end