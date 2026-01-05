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

addpath(genpath('./all_libs/TIP-Code-main_HTNN'))

%%

for iii= 1:6
    a1 = centers(1,1,iii); b1 = centers(1,2,iii);
    a2 = centers(2,1,iii); b2 = centers(2,2,iii);
    a3 = centers(3,1,iii); b3 = centers(3,2,iii);
    a4 = centers(4,1,iii); b4 = centers(4,2,iii);

    wave_data_ = wave_data;
    wave_data_(a1:a4,b1:b2,:) = 0;
    
    Trans = ones(size(wave_data));
    Trans(a1:a4,b1:b2,:) = 0;

    %% Low-Rank High-Order Tensor Completion With Applications in Visual Data  
    %  Wenjin Qin (HTNN-DCT)

    wave_data_ = wave_data;
    wave_data_(a1:a4,b1:b2,:) = 0;
    
    Trans = ones(size(wave_data));
    Trans(a1:a4,b1:b2,:) = 0;

    pos = [];

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
      fprintf('===== t-SVD by Discrete Cosine Transform =====\n');
      U{1}=sqrt(n1)*dct(eye(n1));U{2}=sqrt(n2)*dct(eye(n2));U{3}=sqrt(n3)*dct(eye(n3));
      %U{1}=dct(eye(n3));U{2}=dct(eye(n4));
       t0=tic;
         [Xhat,~,~ ] =HTNN_U(U,XT,mark,opts);
       time = toc(t0);
       Xhat1=max(0,Xhat);
       Xhat2=min(maxP,Xhat1);
       Xhat2=permute(Xhat2,shift_dim);
    
       all_recon_HTNN_DCT{iii} = Xhat2;
        
        wave_power_HTNN_DCT = sum(abs(Xhat2(a1:a4,b1:b2,1:5)).^2,3);
       
        mx = max(wave_power_HTNN_DCT(:));
        [ii,jj] = find(wave_power_HTNN_DCT==mx);
            
        pos = [pos; a1+ii,b1+jj];
        Positions_HTNN_DCT{iii} = pos;
end

rmpath(genpath('./all_libs/TIP-Code-main_HTNN'))
save('./results/HTNN_DCT_results.mat','all_recon_HTNN_DCT',"Positions_HTNN_DCT")