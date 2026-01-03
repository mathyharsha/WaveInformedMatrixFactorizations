function [phi, sample_grid2d] = under_sampl_mtx(dim1, dim2, N)
% Constructing random under-sampling matrix
% dim1 : dimension in vertical (column) driection
% dim2 : dimension in horizontal (row) driection
% N    : The number of sparse sampling points 

% rng(123123);
p = randperm(dim1*dim2)'; 
p_sample = sort(p(1:N));

% Under-sampling matrix
phi = zeros(N, dim1*dim2);    % Sensing matrix

for i = 1:N
    phi(i, p_sample(i)) = 1;
end

sample_grid = zeros(dim1*dim2,1);
sample_grid(p_sample) = 1;

sample_grid2d = reshape(sample_grid, [dim1, dim2]);

% Visualizing the under-sampling matrix
% figure; imagesc(phi)
% figure; imagesc(sample_grid2d); colormap('gray');

end