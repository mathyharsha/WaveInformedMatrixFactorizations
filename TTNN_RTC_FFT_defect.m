clc
clear

data_dam = load('./data/AL0625_1_dam.mat');
data_undam = load('./data/AL0625_1_undam.mat');

rng(42);

%% Preprocess data and mask to have the required shape for further steps

% Wave data with damage present
wave_data_dam_ = data_dam.AL0625_1_dam;
wave_data_dam = reshape(wave_data_dam_,sqrt(size(wave_data_dam_,1)),sqrt(size(wave_data_dam_,1)),size(wave_data_dam_,2));

wave_data_dam = wave_data_dam(:,:,1:100);

% Wave data without damage

wave_data_undam_ = data_undam.AL0625_1_undam;
wave_data_undam = reshape(wave_data_undam_,sqrt(size(wave_data_undam_,1)),sqrt(size(wave_data_undam_,1)),size(wave_data_undam_,2));

wave_data_undam = wave_data_undam(:,:,1:100);


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

addpath(genpath('./all_libs/TTNN-RTC'));

indx = 3;

for indx = 1:6

a1 = centers(1,1,indx); b1 = centers(1,2,indx);
a2 = centers(2,1,indx); b2 = centers(2,2,indx);
a3 = centers(3,1,indx); b3 = centers(3,2,indx);
a4 = centers(4,1,indx); b4 = centers(4,2,indx);


% Damaged

        wave_data_ = wave_data_dam;
        wave_data_(a1:a4,b1:b2,:) = 0;
        
        Trans = ones(size(wave_data_dam));
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
        reconstructedData_dam = X;


% Undamaged

        wave_data_ = wave_data_undam;
        wave_data_(a1:a4,b1:b2,:) = 0;
        
        Trans = ones(size(wave_data_undam));
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
        reconstructedData_undam = X;

%%
      defect = reconstructedData_dam - reconstructedData_undam;
      
      wave_power_est = sum(defect(a1:a4,b1:b2,8:15).^2,3);

      mx = max(wave_power_est(:));
      [ii,jj] = find(wave_power_est==mx);

      TTNN_RTC_FFT_defect_pos = [a1+ii,b1+jj];

      disp_wave_power_est = zeros(size(defect,1),size(defect,2));
      
      disp_wave_power_est(a1:a4,b1:b2) = wave_power_est;

      Positions_defect{indx} = TTNN_RTC_FFT_defect_pos;
      WavePower_est{indx} = disp_wave_power_est;

end
save('./results/TTNN_RTC_FFT_defect.mat','Positions_defect',"WavePower_est")
rmpath(genpath('./all_libs/TTNN-RTC'));