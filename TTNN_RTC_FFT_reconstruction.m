clear all

data = load('./data/AL0625_1_dam.mat');


%% Preprocess data and mask to have the required shape for further steps


wave_data_ = data.AL0625_1_dam;
wave_data = reshape(wave_data_,sqrt(size(wave_data_,1)),sqrt(size(wave_data_,1)),size(wave_data_,2));

wave_power = sum(wave_data.^2,3);
wave_data = wave_data(:,:,1:100);

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

%%
% % Code for the paper: Guangjing Song, Michael K. Ng, and Xiongjun Zhang. 
% % Robust Tensor Completion Using Transformed Tensor Singular Value Decomposition, 
% % Numerical Linear Algebra with Applications, 27(3):e2299, 2020.

addpath(genpath('./all_libs/TTNN-RTC'));

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
        mark = Trans;
        ZZ = wave_data_/max(wave_data_(:));

        dim = size(ZZ);
        %% initial parameters
        beta = 0.05;
        gamma = 1.618;
        MaxIte = 500;
        tol = 5e-4;


        mu = 2/sqrt(max(dim(1),dim(2))*dim(3));
        
        
        %% Initial point
        
        X0 = zeros(dim(1),dim(2),dim(3));
        Y0 = X0;
        Z0 = X0;
        E0 = X0;
        
        
        opts.gamma = gamma;
        opts.tol = tol;   
        
        opts.beta = beta;  
        opts.mu = mu;
        % opts.beta2 = beta2;  
        opts.MaxIte = MaxIte;                   
        opts.X0 = X0;  
        opts.Y0 = Y0;  
        opts.Z0 = Z0;     
        opts.E0 = E0;  
        opts.Omega = mark;          
        opts.dim = dim;    


        %%  FFT
          fprintf('===== t-SVD by Fourier transform =====\n');
           tic
            [X, E, k, eta, aaa] = TNN(ZZ,mark,dim,opts);
           toc
         
        %%

        all_recon_TTNN_RTC{iii} = X;
        
        wave_power_TTNN_RTC = sum(X(a1:a4,b1:b2,1:5).^2,3);
        mx = max(wave_power_TTNN_RTC(:));
        [ii,jj] = find(wave_power_TTNN_RTC==mx);
            
        pos = [pos; a1+ii,b1+jj];
        Positions_TTNN_RTC{iii} = pos;
         
end

rmpath(genpath('./all_libs/TTNN-RTC'));
save('./results/TTNN_RTC_results.mat','all_recon_TTNN_RTC',"Positions_TTNN_RTC")