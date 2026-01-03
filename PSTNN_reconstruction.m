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
pos = [];
%% PSTNN - WORKS!!

for iii= 1:5

    addpath(genpath('./all_libs/Multi-dimensional-imaging-data-recovery-via-minimizing-the-partial-sum-of-tubal-nuclear-norm-master')); 

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
    
    all_recon_PSTNN{iii} = Xhat3;
    
    wave_power_PSTNN = sum(Xhat3(a1:a4,b1:b2,1:5).^2,3);
    mx = max(wave_power_PSTNN(:));
    [ii,jj] = find(wave_power_PSTNN==mx);
        
    pos = [pos; a1+ii,b1+jj];
        
        
    Positions_PSTNN{iii} = pos;

    rmpath(genpath('./all_libs/Multi-dimensional-imaging-data-recovery-via-minimizing-the-partial-sum-of-tubal-nuclear-norm-master')); 
    

end