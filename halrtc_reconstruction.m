clc 
clear all



data = load('./all_libs/AL0625_1_dam.mat');


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
%%

for iii= 1:1
       
    addpath(genpath('./all_libs/TensorCompletion'))
    
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
    all_recons_HaLRTC{iii} = reconstructedData_HaLRTC;
    wave_power_HaLRTC = sum(reconstructedData_HaLRTC(a1:a4,b1:b2,1:5).^2,3);
    mx = max(wave_power_HaLRTC(:));
    [ii,jj] = find(wave_power_HaLRTC==mx);
    
    pos = [pos; a1+ii,b1+jj];
    
    
    Positions_HaLRTC{iii} = pos;

    rmpath(genpath('./all_libs/TensorCompletion'))


end
