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

    %% Low-Rank High-Order Tensor Completion With Applications in Visual Data  
    %  Wenjin Qin (HTNN-DCT)
    addpath(genpath('./all_libs/TIP-Code-main_HTNN'))

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

    

    X = wave_data_/max(wave_data_(:));
    mark = Trans;
    maxP = max(abs(wave_data_(:)));
    
    [d1, d2 ,d3] = size(X);
    
    shift_dim=[1,2,3];
    XT = permute(X,shift_dim);
    mark = permute(mark,shift_dim);
    [n1, n2 ,n3] = size(XT);
    
%% initial parameters 
opts.mu = 1e-4;
opts.max_mu = 1e8;
opts.max_iter =80;
opts.DEBUG = 1;
opts.rho = 1.2;
opts.tol = 1e-10;
%opts.tol = 1e-8;
    
    %%  Discrete Cosine Transform (DCT)
       fprintf('===== t-SVD by FFT =====\n');
    t0=tic;
     [Xhat,~,~ ] =HTNN_FFT(XT,mark,opts);
   time = toc(t0);
       Xhat1=max(0,Xhat);
       Xhat2=min(maxP,Xhat1);
       Xhat2=permute(Xhat2,shift_dim);
    
       all_recon_HTNN_FFT{iii} = Xhat2;
        
        wave_power_HTNN_FFT = sum(abs(Xhat2(a1:a4,b1:b2,1:5)).^2,3);
       
        mx = max(wave_power_HTNN_FFT(:));
        [ii,jj] = find(wave_power_HTNN_FFT==mx);
            
        pos = [pos; a1+ii,b1+jj];
            
            
        Positions_HTNN_FFT{iii} = pos;
        
  
       rmpath(genpath('./all_libs/TIP-Code-main_HTNN'))


end

