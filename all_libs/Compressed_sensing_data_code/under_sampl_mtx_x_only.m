function [phi, sample_grid2d] = under_sampl_mtx_x_only(dim1, dim2, comp_rate)

% number of points in a horizontal line scan
N_und_line = round(dim2*(1-comp_rate));

p = randperm(dim2); p_sample = sort(p(1:N_und_line))';

sample_grid2d = zeros(dim1, dim2);

sample_grid2d(:, p_sample) = 1;

% Convret grid 2D matrix to down-sampling matrix
phi = grid2d_to_phi(sample_grid2d);
end