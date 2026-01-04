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

        addpath(genpath('./all_libs/TTNN-RTC'));
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
           %% U Transform

        fprintf('===== t-SVD by unitary transform =====\n');
        
        %% U matrix
        O = tenmat(X,[3]); % unfolding
        O = O.data;
        [U, D, V] = svd(O,'econ');
        
        
        mu = 25/sqrt(max(dim(1),dim(2))*dim(3));
        opts.mu = mu;
        %%
          tic
            [XU, E, k] = UTNN(U, ZZ,mark,dim,opts);
          toc
         
        %%

        all_recon_TTNN_RTC_U{iii} = XU;
        
        wave_power_TTNN_RTC_U = sum(XU(a1:a4,b1:b2,1:5).^2,3);
        mx = max(wave_power_TTNN_RTC_U(:));
        [ii,jj] = find(wave_power_TTNN_RTC_U==mx);
            
        pos = [pos; a1+ii,b1+jj];
            
            
        Positions_TTNN_RTC_U{iii} = pos;
         
        rmpath(genpath('./all_libs/TTNN-RTC'));


end