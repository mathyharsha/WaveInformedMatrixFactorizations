





addpath("C:\Users\SmartDATA Lab\Documents\MATLAB\ICCV23\WIMF\WIMF\")


clc

data = load('AL0625_1_dam.mat');


%% Preprocess data as required for further steps and compute Laplacian and its eigenvalues and eigenvectors

wave_data_str = data.AL0625_1_dam;
wave_data = reshape(wave_data_str,sqrt(size(wave_data_str,1)),sqrt(size(wave_data_str,1)),size(wave_data_str,2));

wave_data = wave_data(:,:,1:50); 

wave_power = sum(wave_data.^2,3);

mx = max(wave_power(:));
[trueI,trueJ] = find(wave_power==mx);

wave_data_ = wave_data;
wave_data_(6:90,20:90,:)=0;


X = wave_data_;


Trans = ones(size(wave_data_));
Trans(6:90,20:90,:) = 0;


omega = find(Trans==1);


opts.mu = 10^-3;
opts.tol = 1e-7;
opts.rho = 1.5;
opts.max_iter = 200;
opts.DEBUG = 0;
opts.max_mu = 1e10;

PSNR0 = psnr(wave_data_,wave_data);
SSIM0 = ssim(wave_data_,wave_data);

%% Multi-dimensional imaging data recovery via minimizing the partial sum of tubal nuclear norm 
% Journal of Computational and Applied Mathematics 2020, Tai-Xiang Jiang et al.


rho = 1;% tune this parameter to controal the estimated rank N
[rankN,~] = prox_rankN(wave_data,rho);%n*ones(1,n3);%
[Xhat3,~,~,~] =  LRTC_PSTNN(wave_data_,omega,opts,rankN);%,Xhat2);%

Xhat3 = max(Xhat3,0);
Xhat3 = min(Xhat3,maxP);
for i = 1:size(wave_data,3)
SSIMv3(i) = ssim(Xhat3(:,:,i),X(:,:,i));
PSNRv3(i) = psnr(Xhat3(:,:,i),X(:,:,i));
end
mSSIM3 = mean(SSIMv3);mPSNR3 = mean(PSNRv3);