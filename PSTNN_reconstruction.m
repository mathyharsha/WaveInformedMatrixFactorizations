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


addpath(genpath('./all_libs/Multi-dimensional-imaging-data-recovery-via-minimizing-the-partial-sum-of-tubal-nuclear-norm-master')); 
%% PSTNN - WORKS!!
% Tai-Xiang Jiang, Ting-Zhu Huang, Xi-Le Zhao, and Liang-Jian Deng. 
% Multi-dimensional imaging data recovery via minimizing the partial
% sum of tubal nuclear norm.

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

    
    omega = find(Trans==1);
    
    
    opts.mu = 10^-3;
    opts.tol = 1e-7;
    opts.rho = 1.5;
    opts.max_iter = 200;
    opts.DEBUG = 0;
    opts.max_mu = 1e10;
    
    
    maxP = max(abs(wave_data(:)));
    
    rho = 1;% tune this parameter to control the estimated rank N
    [rankN,~] = prox_rankN(wave_data,rho);%n*ones(1,n3);%
    [Xhat3,~,~,~] =  LRTC_PSTNN(wave_data_,omega,opts,rankN);%,Xhat2);%
    
    Xhat3 = max(Xhat3,0);
    Xhat3 = min(Xhat3,maxP);
    % for i = 1:size(wave_data,3)
    % SSIMv3(i) = ssim(Xhat3(:,:,i),X(:,:,i));
    % PSNRv3(i) = psnr(Xhat3(:,:,i),X(:,:,i));
    % end
    % mSSIM3 = mean(SSIMv3);mPSNR3 = mean(PSNRv3);
    
    all_recons_PSTNN{iii} = Xhat3;
    
    wave_power_PSTNN = sum(Xhat3(a1:a4,b1:b2,1:5).^2,3);
    mx = max(wave_power_PSTNN(:));
    [ii,jj] = find(wave_power_PSTNN==mx);
        
    pos = [pos; a1+ii,b1+jj];
    Positions_PSTNN{iii} = pos;

end

rmpath(genpath('./all_libs/Multi-dimensional-imaging-data-recovery-via-minimizing-the-partial-sum-of-tubal-nuclear-norm-master')); 
save('./results/PSTNN_results.mat','all_recons_PSTNN',"Positions_PSTNN")