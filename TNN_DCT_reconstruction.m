%% Low-Rank Tensor Completion With a New Tensor Nuclear Norm Induced by Invertible
% Linear Transforms IEEE Conference on Computer Vision and Pattern Recognition 2019,
% Canyi Lu et al. [TNN-DCT] -- WORKS!

clc 
clear all



data = load('./data/AL0625_1_dam.mat');


%% Preprocess data and mask to have the required shape for further steps

wave_data_ = data.AL0625_1_dam;
wave_data = reshape(wave_data_,sqrt(size(wave_data_,1)),sqrt(size(wave_data_,1)),size(wave_data_,2));

wave_data = wave_data(:,:,1:100);



wave_power = sum(wave_data.^2,3);

mx = max(wave_power(:));
[trueI,trueJ] = find(wave_power==mx);

% Lx = 4*del2(eye(size(wave_data,2)));
% Ly = 4*del2(eye(size(wave_data,1)));
% Lt = 4*del2(eye(size(wave_data,3)));
% 
% [V_x,Dx] = eig(Lx);
% [V_y,Dy] = eig(Ly);
% [V_t,Dt] = eig(Lt);

centers = [1;1;1;1]*[trueI,trueJ];
pos = [];
%%

for iii= 1:5
       
    addpath(genpath('./all_libs/Tensor-robust-PCA-and-tensor-completion-under-linear-transform-main'))
    
    [a1,b1] = mids(1,1,centers(1,1),centers(1,2));
    [a2,b2] = mids(1,100,centers(2,1),centers(2,2));
    [a3,b3] = mids(100,100,centers(3,1),centers(3,2));
    [a4,b4] = mids(100,1,centers(4,1),centers(4,2));
    
    centers(1,1) = a1; centers(1,2) = b1;
    centers(2,1) = a2; centers(2,2) = b2;
    centers(3,1) = a3; centers(3,2) = b3;
    centers(4,1) = a4; centers(4,2) = b4;
    wave_data_ = wave_data;
    wave_data_(a1:a4,b1:b2,:) = 0;
    
    Trans = ones(size(wave_data));
    Trans(a1:a4,b1:b2,:) = 0;
    [n1,n2,n3] = size(wave_data_);
    r = 10;
    % maxP = max(abs(wave_data_(:)));
    % dimX = size(wave_data_);
    % 
    % M2 = Frontal2Lateral(wave_data_); % each lateral slice is a channel of the image
    % omega2 = zeros(dimX);
    % Iones = ones(dimX);
    % omega2 = Trans;
    % omega2 = Frontal2Lateral(omega2);
    % omega2 = find(omega2==1);
    % n3 = size(M2,3);
    % 
    % % transform.L = @fft; transform.l = n3; transform.inverseL = @ifft;
    % transform.L = @dct; transform.l = 1; transform.inverseL = @idct;
    % % L = dftmtx(n3); transform.l = n3; transform.L = L;
    % % L = dct(eye(n3)); transform.l = 1; transform.L = L;
    % % L = RandOrthMat(n3); transform.l = 1; transform.L = L;
    % 
    % opts.DEBUG = 1;
    % Xhat2 = lrtc_tnn(M2,omega2,transform,opts);
    % % Xhat2 = max(Xhat2,0);
    % % Xhat2 = min(Xhat2,maxP);
    % Xhat2 = Lateral2Frontal(Xhat2);
    % reconstructedData_TNN_DCT = Xhat2;
    % all_recons_TNN_DCT{iii} = reconstructedData_TNN_DCT;
    % wave_power_TNN_DCT = sum(reconstructedData_TNN_DCT(a1:a4,b1:b2,1:5).^2,3);
    % mx = max(wave_power_TNN_DCT(:));
    % [ii,jj] = find(wave_power_TNN_DCT==mx);
    % 
    % pos = [pos; a1+ii,b1+jj];
    % 
    % 
    % Positions_TNN_DCT{iii} = pos;

    L = RandOrthMat(n3); transform.l = 1; transform.L = @dct; transform.inverseL = @idct;

    X = wave_data_;
    trankX = tubalrank(X,transform);
    
    dr = (n1+n2-r)*r*n3;
    rho = 3;
    m = rho*dr;
    p = m/(n1*n2*n3);
    
    omega = find(Trans==1);
    M = zeros(n1,n2,n3);
    M(omega) = X(omega);
    
    opts.DEBUG = 1;
    Xhat = lrtc_tnn(M,omega,transform,opts);
    
    trank = tubalrank(Xhat,transform);
    RSE = norm(X(:)-Xhat(:))/norm(X(:));
    
    all_recon_TNN_DCT{iii} = Xhat;
    
    wave_power_TNN_DCT = sum(Xhat(a1:a4,b1:b2,1:5).^2,3);
    mx = max(wave_power_TNN_DCT(:));
    [ii,jj] = find(wave_power_TNN_DCT==mx);
        
    pos = [pos; a1+ii,b1+jj];
        
        
    Positions_TNN_DCT{iii} = pos;


    rmpath(genpath('./all_libs/Tensor-robust-PCA-and-tensor-completion-under-linear-transform-main'))

end
