clc
clear all;

% addpath('C:\Users\SmartDATALab\Documents\MATLAB\ICCV_23\Compressed_sensing_data_code');
% addpath('C:\Users\SmartDATALab\Documents\MATLAB\ICCV_23\spgl1-2.1');
% addpath('C:\Users\SmartDATALab\Documents\MATLAB\ICCV_23\TensorCompletion\mylib');
% addpath('C:\Users\SmartDATALab\Documents\MATLAB\ICCV_23\TensorCompletion')
% addpath('C:\Users\SmartDATALab\Documents\MATLAB\ICCV_23\Tensor_Completion_and_Tensor_RPCA-master')
% addpath('C:\Users\SmartDATALab\Documents\MATLAB\ICCV_23\Tensor_Completion_and_Tensor_RPCA-master\tSVD')
% addpath('C:\Users\SmartDATALab\Documents\MATLAB\ICCV_23\Tensor_Completion_and_Tensor_RPCA-master\solvers')
% addpath('C:\Users\SmartDATALab\Documents\MATLAB\ICCV_23\Tensor_Completion_and_Tensor_RPCA-master\proxFunctions')



data = load('AL0625_1_dam.mat');


%% Preprocess data and mask to have the required shape for further steps

wave_data_ = data.AL0625_1_dam;
wave_data = reshape(wave_data_,sqrt(size(wave_data_,1)),sqrt(size(wave_data_,1)),size(wave_data_,2));

wave_data = wave_data(:,:,1:100);



wave_power = sum(wave_data.^2,3);

mx = max(wave_power(:));
[trueI,trueJ] = find(wave_power==mx);

Lx = 4*del2(eye(size(wave_data,2)));
Ly = 4*del2(eye(size(wave_data,1)));
Lt = 4*del2(eye(size(wave_data,3)));

[V_x,Dx] = eig(Lx);
[V_y,Dy] = eig(Ly);
[V_t,Dt] = eig(Lt);

centers = [1;1;1;1]*[trueI,trueJ];
%%

for iii= 1:5

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

%% WIMF
[Dc,x,C] = waveInformedMatFac(wave_data_,'count',5);

all_Ds_WIMFs{iii} = Dc;
pos = [];
for jjj=1:size(Dc,2)
    
    power_D_WIMF = inv_spectral_wave_transform(Dc(:,jjj),V_x,V_y,V_t);
    wave_power_WIMF = sum(power_D_WIMF(a1:a4,b1:b2,1:5).^2,3);
    mx = max(wave_power_WIMF(:));
    [ii,jj] = find(wave_power_WIMF==mx);
    pos = [pos; a1+ii,b1+jj];
end


reconstructedData_WIMF = inv_spectral_wave_transform(Dc*x,V_x,V_y,V_t);
all_recons_WIMF{iii} = reconstructedData_WIMF;
wave_power_WIMF = sum(reconstructedData_WIMF(a1:a4,b1:b2,1:5).^2,3);
mx = max(wave_power_WIMF(:));
[ii,jj] = find(wave_power_WIMF==mx);



pos = [pos; a1+ii,b1+jj];
    
    
Positions_WIMF{iii} = pos;

%% SPG-MMV with DCT
% Development of an MEMS ultrasonic microphone array system and its
% application to compressed wavefield imaging of concrete, Homin Song et al.

pos = [];
WF_original = wave_data_;

dim1 = size(WF_original, 1) ; dim2 = size(WF_original, 2); dim3 = size(WF_original, 3);   
% 
m = dim1*dim2;          % total number of true data within the grid
comp_rate = 0.85;    % Compression rate
k = round(m*(1-comp_rate));                % the number of measurements within a pixlated grid

% Under-sampling matrix: phi
% [phi, sample_grid2d] = under_sampl_mtx(dim1, dim2, k);

%[phi, sample_grid2d] = under_sampl_mtx_x_only(dim1, dim2, comp_rate);

sample_grid2d = ones(dim1,dim2);
sample_grid2d(a1:a4,b1:b2,:) = 0;
phi = grid2d_to_phi(sample_grid2d);

% Basis matrix: B
B = dct2_basis_mtx(dim1, dim2);

% Sensing matrix: Psi
Psi = phi*B;

%Input for SPG-MMV (DCT)
Yt = phi*reshape(WF_original, [dim1*dim2, dim3]);   % Sparse measurements
n = 40;                  % n should be smaller than the size of WF_original in dim3
% Yf = dct(Yt, n, 2);     % dct(WF_original, n, dim)
Yf = dct(Yt, dim3, 2);
Yf = Yf(:,1:n);

% Solving MMV problem using SPG-MMV
opts = spgSetParms('verbosity',0);         % Turn off the SPGL1 log output
A = spg_mmv(Psi, Yf, 0, opts);             % Reconstructed basis coefficient matrix A
Xf = B*A;                                  % Reconstruction in frequency domain

% Converting into time domain
Xt = idct(Xf, dim3, 2); 

% Reshaping the reconstructed wavefield
Xt = reshape(Xt, [dim1, dim2, dim3]);

psnr = 20*log10(sqrt(dim1*dim2*dim3)*max(abs(WF_original(:)))/norm(reshape(WF_original-Xt, [dim1*dim2*dim3, 1]),2));
disp(psnr)    


reconstructedData_MMV_1 = Xt;
all_recons_MMV_1{iii} = reconstructedData_MMV_1;
wave_power_MMV_1 = sum(reconstructedData_MMV_1(a1:a4,b1:b2,1:5).^2,3);
mx = max(wave_power_MMV_1(:));
[ii,jj] = find(wave_power_MMV_1==mx);

pos = [pos; a1+ii,b1+jj];
    
Positions_MMV_1{iii} = pos;



%% SPG-MMV without DCT
% Accelerated non-contact guided wave array imaging via sparse array data
% reconstruction
pos = [];
WF_original = wave_data_;

dim1 = size(WF_original, 1) ; dim2 = size(WF_original, 2); dim3 = size(WF_original, 3);   
% 
m = dim1*dim2;          % total number of true data within the grid
comp_rate = 0.85;    % Compression rate
k = round(m*(1-comp_rate));                % the number of measurements within a pixlated grid

% Under-sampling matrix: phi
% [phi, sample_grid2d] = under_sampl_mtx(dim1, dim2, k);

%[phi, sample_grid2d] = under_sampl_mtx_x_only(dim1, dim2, comp_rate);

sample_grid2d = ones(dim1,dim2);
%sample_grid2d(36:75,36:75,:) = 0; 
sample_grid2d(6:90,20:90,:) = 0;
phi = grid2d_to_phi(sample_grid2d);

% Basis matrix: B
B = dct2_basis_mtx(dim1, dim2);

% Sensing matrix: Psi
Psi = phi*B;

%Input for SPG-MMV (DCT)
Yt = phi*reshape(WF_original, [dim1*dim2, dim3]);   % Sparse measurements
n = dim3;                  % n should be smaller than the size of WF_original in dim3
% Yf = dct(Yt, n, 2);     % dct(WF_original, n, dim)
% Yf = dct(Yt, dim3, 2);
% Yf = Yf(:,1:n);

% Solving MMV problem using SPG-MMV
opts = spgSetParms('verbosity',0);         % Turn off the SPGL1 log output
A = spg_mmv(Psi, Yt, 0, opts);             % Reconstructed basis coefficient matrix A
Xt = B*A;                                  % Reconstruction in frequency domain

% Converting into time domain
%Xt = idct(Xf, dim3, 2); 

% Reshaping the reconstructed wavefield
Xt = reshape(Xt, [dim1, dim2, dim3]);

psnr = 20*log10(sqrt(dim1*dim2*dim3)*max(abs(WF_original(:)))/norm(reshape(WF_original-Xt, [dim1*dim2*dim3, 1]),2));

reconstructedData_MMV_2 = Xt;
all_recons_MMV_2{iii} = reconstructedData_MMV_2;
wave_power_MMV_2 = sum(reconstructedData_MMV_2(a1:a4,b1:b2,1:5).^2,3);
mx = max(wave_power_MMV_2(:));
[ii,jj] = find(wave_power_MMV_2==mx);

pos = [pos; a1+ii,b1+jj];
    
Positions_MMV_2{iii} = pos;
    

trans = reshape(Trans,size(Trans,1)*size(Trans,2)*size(Trans,3),1);

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

%% Another tensor completion
% Novel methods for multilinear data completion and de-noising based on
% tensor-SVD, Zemin Zhang et. al.

pos = [];
trans = reshape(Trans,size(Trans,1)*size(Trans,2)*size(Trans,3),1);

normalize = max(wave_data(:));

b = wave_data(:)/max(wave_data(:));


[n1,n2,n3]             =        size(wave_data)               ;
alpha                  =        1                             ;
maxItr                 =        1000                          ; % maximum iteration
rho                    =        0.01                          ;
myNorm                 =        'tSVD_1'                      ; % dont change for now

%================ main process of completion =======================
X   =    tensor_cpl_admm( trans , b , rho , alpha , ...
                     [n1,n2,n3] , maxItr , myNorm , 0 );
X                      =        X * normalize                 ;
X                      =        reshape(X,[n1,n2,n3])         ;
            
X_dif                  =        wave_data-X                           ;
RSE                    =        norm(X_dif(:))/norm( wave_data(:))     ;

reconstructedData_TC  =X;
all_recons_TC{iii} = reconstructedData_TC;
wave_power_TC = sum(reconstructedData_TC(a1:a4,b1:b2,1:5).^2,3);
mx = max(wave_power_TC(:));
[ii,jj] = find(wave_power_TC==mx);

pos = [pos; a1+ii,b1+jj];
    
Positions_TC{iii} = pos;





end

