function basis_vec = dct2_basis_mtx(dim1, dim2)
% Vectorized 2-D DCT basis matrix construction

B1 = dctmtx(dim1)';
B2 = dctmtx(dim2)';

basis = zeros(dim1,dim2,dim1*dim2);
basis_num = 1;

for i=1:dim1
    for j=1:dim2
            basis(:,:,basis_num) = B1(:,i)*B2(:,j)';
            basis_num = basis_num + 1;
    end
end

basis_vec = reshape(basis, [dim1*dim2, dim1*dim2]);
end