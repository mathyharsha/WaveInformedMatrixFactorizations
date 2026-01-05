%% Low-Rank Tensor Completion With a New Tensor Nuclear Norm Induced by Invertible
% Linear Transforms IEEE Conference on Computer Vision and Pattern Recognition 2019,
% Canyi Lu et al. [TNN-DCT] -- WORKS!

clc 
clear all



data = load('./data/AL0625_1_dam.mat');


%% Preprocess data and mask to have the required shape for further steps

wave_data_ = data.AL0625_1_dam;
wave_data = reshape(wave_data_,sqrt(size(wave_data_,1)),sqrt(size(wave_data_,1)),size(wave_data_,2));

wave_power = sum(wave_data.^2,3);

wave_data = wave_data(:,:,1:100);

centers = zeros(4,2,6);

A1B1 = floor([0.15,0.20;0.06,0.20;0.10,0.10;0.06,0.06;0.04,0.04;0.02,0.02]*100);
A3B3 = floor([0.83,0.85;0.90,0.90;0.90,0.90;0.94,0.94;0.96,0.96;0.98,0.98]*100);

for i = 1:6
    a1 = A1B1(i,1);
    b1 = A1B1(i,2);
    a3 = A3B3(i,1);
    b3 = A3B3(i,2);
    
    a2 = a1;
    b2 = b3;
    b4 = b1;
    a4 = a3;

    centers(1,1,i) = a1; centers(1,2,i) = b1;
    centers(2,1,i) = a2; centers(2,2,i) = b2;
    centers(3,1,i) = a3; centers(3,2,i) = b3;
    centers(4,1,i) = a4; centers(4,2,i) = b4;
end


mx = max(wave_power(:));
[trueI,trueJ] = find(wave_power==mx);

addpath(genpath('./all_libs/Tensor-robust-PCA-and-tensor-completion-under-linear-transform-main'))
    

%%

for iii= 1:6
       
    pos = [];
    
    a1 = centers(1,1,iii); b1 = centers(1,2,iii);
    a2 = centers(2,1,iii); b2 = centers(2,2,iii);
    a3 = centers(3,1,iii); b3 = centers(3,2,iii);
    a4 = centers(4,1,iii); b4 = centers(4,2,iii);

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
end

rmpath(genpath('./all_libs/Tensor-robust-PCA-and-tensor-completion-under-linear-transform-main'))
save('./results/TNN_DCT_results.mat','all_recon_TNN_DCT',"Positions_TNN_DCT")