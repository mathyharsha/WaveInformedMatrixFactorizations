clc 
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


addpath(genpath('./all_libs/TensorCompletion'))

for iii= 1:6

    a1 = centers(1,1,iii); b1 = centers(1,2,iii);
    a2 = centers(2,1,iii); b2 = centers(2,2,iii);
    a3 = centers(3,1,iii); b3 = centers(3,2,iii);
    a4 = centers(4,1,iii); b4 = centers(4,2,iii);

    wave_data_ = wave_data;
    wave_data_(a1:a4,b1:b2,:) = 0;
    
    Trans = ones(size(wave_data));
    Trans(a1:a4,b1:b2,:) = 0;

    %% HaLRTC Tensor Completion
    % Tensor Completion for Estimating Missing Values in Visual Data
    % Ji Liu et al.
    
    T = wave_data;
    Omega = logical(Trans);
    
    alpha = [10, 10, 1e-3];
    alpha = alpha / sum(alpha);
    
    maxIter = 500;
    epsilon = 1e-5;
    
    % "X" returns the estimation, 
    % "errList" returns the list of the relative difference of outputs of two neighbor iterations 
    
    % High Accuracy LRTC (solve the original problem, HaLRTC algorithm in the paper)
    
    pos = [];
    pos_secondary = [];
    
    rho = 1e-1;
    [X_H, errList_H] = HaLRTC(...
        T, ...                       % a tensor whose elements in Omega are used for estimating missing value
        Omega,...               % the index set indicating the obeserved elements
        alpha,...                  % the coefficient of the objective function, i.e., \|X\|_* := \alpha_i \|X_{i(i)}\|_* 
        rho,...                      % the initial value of the parameter; it should be small enough  
        maxIter,...               % the maximum iterations
        epsilon...                 % the tolerance of the relative difference of outputs of two neighbor iterations 
        );
    
    reconstructedData_HaLRTC  =X_H;
    psnr_vals{iii} = calculate_psnr(wave_data,reconstructedData_HaLRTC);
    all_recons{iii} = reconstructedData_HaLRTC;
    wave_power_HaLRTC = sum(reconstructedData_HaLRTC(a1:a4,b1:b2,1:5).^2,3);
    mx = max(wave_power_HaLRTC(:));
    [ii,jj] = find(wave_power_HaLRTC==mx);
    
     mg = min(wave_power(:));
     [ij,ji] = find(wave_power==mg);
        
    
    pos = [pos; a1+ii,b1+jj];
    pos_secondary = [pos_secondary; a1+ij,b1+ji];
    Positions{iii} = pos;
    Positions_secondary{iii} = pos_secondary;


end

rmpath(genpath('./all_libs/TensorCompletion'))
save('./results/HALRTC_results.mat','all_recons',"Positions","Positions_secondary","psnr_vals")

