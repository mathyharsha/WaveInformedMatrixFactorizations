clc

% addpath('C:\Users\SmartDATALab\Documents\MATLAB\ICCV_23\TIP-Code-main_HTNN\'); 
% addpath('C:\Users\SmartDATALab\Documents\MATLAB\ICCV_23\Tensor-robust-PCA-and-tensor-completion-under-linear-transform-main\');
% addpath('C:\Users\SmartDATALab\Documents\MATLAB\ICCV_23\Tensor_Completion_and_Tensor_RPCA-master\')
% addpath('C:\Users\SmartDATALab\Documents\MATLAB\ICCV_23\Tensor_Completion_and_Tensor_RPCA-master\tSVD')
% addpath('C:\Users\SmartDATALab\Documents\MATLAB\ICCV_23\Tensor_Completion_and_Tensor_RPCA-master\solvers')
% addpath('C:\Users\SmartDATALab\Documents\MATLAB\ICCV_23\Tensor_Completion_and_Tensor_RPCA-master\pr')
% addpath('C:\Users\SmartDATALab\Documents\MATLAB\ICCV_23\Multi-dimensional-imaging-data-recovery-via-minimizing-the-partial-sum-of-tubal-nuclear-norm-master\'); 
% 

addpath(genpath('/Users/harshat/Downloads/UF_work/all_libs/TIP-Code-main_HTNN'))
addpath(genpath('/Users/harshat/Downloads/UF_work/all_libs/Tensor_Completion_and_Tensor_RPCA-master'))
addpath(genpath('/Users/harshat/Downloads/UF_work/all_libs/Tensor-robust-PCA-and-tensor-completion-under-linear-transform-main'))
addpath(genpath('/Users/harshat/Downloads/UF_work/all_libs/Multi-dimensional-imaging-data-recovery-via-minimizing-the-partial-sum-of-tubal-nuclear-norm-master')); 


data = load('AL0625_1_dam.mat');




%% Preprocess data as required for further steps and compute Laplacian and its eigenvalues and eigenvectors

wave_data_str = data.AL0625_1_dam;
wave_data = reshape(wave_data_str,sqrt(size(wave_data_str,1)),sqrt(size(wave_data_str,1)),size(wave_data_str,2));

wave_data = wave_data(:,:,1:300); 

wave_power = sum(wave_data.^2,3);

mx = max(wave_power(:));
[trueI,trueJ] = find(wave_power==mx);

centers = [1;1;1;1]*[trueI,trueJ];

Lx = 4*del2(eye(size(wave_data,2)));
Ly = 4*del2(eye(size(wave_data,1)));
Lt = 4*del2(eye(size(wave_data,3)));

[V_x,Dx] = eig(Lx);
[V_y,Dy] = eig(Ly);
[V_t,Dt] = eig(Lt);
%%

for iii= 1:1

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

% 
% [Dc,x,C] = waveInformedMatFac(wave_data_);
% 
% 
% all_Ds{iii} = Dc;
% pos = [];
% for jjj=1:size(Dc,2)
%     
%     power_D = inv_spectral_wave_transform(Dc(:,jjj),V_x,V_y,V_t);
%     wave_power = sum(power_D(a1:a4,b1:b2,1:5).^2,3);
%     mx = max(wave_power(:));
%     [ii,jj] = find(wave_power==mx);
%     pos = [pos; a1+ii,b1+jj];
% end
% 
% reconstructedData = inv_spectral_wave_transform(Dc*x,V_x,V_y,V_t);
% wave_power = sum(reconstructedData(a1:a4,b1:b2,1:5).^2,3);
% mx = max(wave_power(:));
% [ii,jj] = find(wave_power==mx);
%     
% pos = [pos; a1+ii,b1+jj];
%     
%     
% Positions{iii} = pos;


%% Multi-dimensional imaging data recovery via minimizing the partial sum of tubal nuclear norm 
% Journal of Computational and Applied Mathematics 2020, Tai-Xiang Jiang et
% al. [PSTNN]

pos = [];
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



%% Low-Rank Tensor Completion With a New Tensor Nuclear Norm Induced by Invertible
% Linear Transforms IEEE Conference on Computer Vision and Pattern Recognition 2019,
% Canyi Lu et al. [TNN-DCT]

opts=  [];

pos = [];

[n1,n2,n3] = size(wave_data_);
r =10; % tubal rank

% transform.L = @fft; transform.l = n3; transform.inverseL = @ifft;
% transform.L = @dct; transform.l = 1; transform.inverseL = @idct;
% L = dftmtx(n3); transform.l = n3; transform.L = L;
% L = dct(eye(n3)); transform.l = 1; transform.L = L;
L = RandOrthMat(n3); transform.l = 1; transform.L = @dct; transform.inverseL = @idct;

X = wave_data_;
trankX = tubalrank(X,transform);

dr = (n1+n2-r)*r*n3;
rho = 3;
m = rho*dr;
p = m/(n1*n2*n3);

omega = find(Trans==1);
M = zeros(n1,n2,n3);
M(omega) = X(omega);

opts.DEBUG = 1;
Xhat = lrtc_tnn(M,omega,transform,opts);

trank = tubalrank(Xhat,transform);
RSE = norm(X(:)-Xhat(:))/norm(X(:));

all_recon_lin_tubal{iii} = Xhat;

wave_power_lin_tubal = sum(Xhat(a1:a4,b1:b2,1:5).^2,3);
mx = max(wave_power_lin_tubal(:));
[ii,jj] = find(wave_power_lin_tubal==mx);
    
pos = [pos; a1+ii,b1+jj];
    
    
Positions_lin_tubal{iii} = pos;


% % %% Framelet Representation of Tensor Nuclear Norm for Third-Order Tensor Completion 
% % % IEEE Transactions on Image Processing 2020, Tai-Xiang Jiang et al. (Not worth showing)
% % 
% %     X = wave_data_;
% %     [n1,n2,n3] = size(X);
% %     for i=1:n3
% %         temp=X(:,:,i);
% %         X(:,:,i)=temp/max(temp(:));
% %     end
% %     
% %     %% Sampling with random position
% %     sample_ratio = 0.3;
% %     fprintf('\n');
% %     fprintf('===========data no. %d=====Results=p=%f======================\n',data_num,sample_ratio);
% %     
% %     Y_tensorT   = X;
% %     Ndim        = ndims(Y_tensorT);
% %     Nway        = size(Y_tensorT);
% %     Omega       = find(Trans==1);
% %     Ind         = zeros(size(X));
% %     Ind(Omega)  = 1;
% %     Y_tensor0   = zeros(Nway);
% %     Y_tensor0(Omega) = Y_tensorT(Omega);
% %     %%
% %     i  = 1;
% %     Re_tensor{i} = Y_tensor0;
% %     [MPSNRALL(i), SSIMALL(i), FSIMALL(i)] = quality(Y_tensorT*255, Re_tensor{i}*255);
% %     time(i) = 0;
% %     enList = 1;
% %     fprintf(' %8.8s    %5.4s    %5.4s    %5.4s   %5.4s   %5.4s \n','method','PSNR', 'SSIM', 'FSIM','iter','time');
% %     fprintf(' %8.8s    %5.3f    %5.3f    %5.3f    %3.3d     %.1f \n',...
% %         methodname{enList(i)},MPSNRALL(enList(i)), SSIMALL(enList(i)), FSIMALL(enList(i)),0,time(i));
% %     %% Perform  algorithm
% %     
% %     
% %     %% Use FTNN
% %     i = i+1;
% %     if EN_FTNN
% %         enList = [enList,i];
% %         %parameters
% %         opts.Frame    = 1; % (0,1,3)
% %         opts.Level    = 4;  % [1,2,3,4,5,6]
% %         opts.wLevel   = -1;
% %         opts.lambda1  = 1;
% %         opts.beta     = 1;
% %         opts.tol      = 1e-2;
% %         opts.rho      = 1;
% %         opts.DEBUG    = 0;
% %         opts.max_iter = 200;
% %         opts.max_beta = 1e10;
% %         
% %         tStart = tic;
% %         [Xhat3_out,obj,~,iter] = LRTC_FL(Y_tensor0,Omega,opts,X);
% %         Re_tensor{i} = Xhat3_out;
% %         time(i)= toc(tStart);
% %         [MPSNRALL(i), SSIMALL(i), FSIMALL(i)] = quality(Y_tensorT*255, Re_tensor{i}*255);
% %         fprintf(' %8.8s    %5.3f    %5.3f    %5.3f    %3.3d   %.1f | Frame =  %d, Level = %d, beta = %.2f\n',...
% %             methodname{i},MPSNRALL(i), SSIMALL(i), FSIMALL(i),iter,time(i),...
% %             opts.Frame,opts.Level,opts.beta);
% %         
% %         
% %     end

%% Low-Rank High-Order Tensor Completion With Applications in Visual Data  
%  Wenjin Qin (HTNN-DCT)

pos = [];

X = wave_data_/max(wave_data_(:));
mark = Trans;


[d1, d2 ,d3,d4] = size(X);

shift_dim=[3,4,1,2];
XT = permute(X,shift_dim);
mark = permute(mark,shift_dim);
[n1, n2 ,n3,n4] = size(XT);

%% initial parameters 
opts.mu = 1e-4;
opts.max_mu = 1e8;
opts.max_iter =500;
opts.DEBUG = 1;
opts.rho = 10;
opts.tol = 1e-10;


%%  Discrete Cosine Transform (DCT)
  fprintf('===== t-SVD by Discrete Cosine Transform =====\n');
  U{1}=sqrt(n3)*dct(eye(n3));U{2}=sqrt(n4)*dct(eye(n4));
 % U{1}=dct(eye(n3));U{2}=dct(eye(n4));
   t0=tic;
     [Xhat,~,~ ] =HTNN_U(U,XT,mark,opts);
   time = toc(t0);
   Xhat1=max(0,Xhat);
   Xhat2=min(maxP,Xhat1);
   Xhat2=permute(Xhat2,shift_dim);

   all_recon_low_high_TIP{iii} = Xhat2;

wave_power_low_high_TIP = sum(Xhat2(a1:a4,b1:b2,1:5).^2,3);
mx = max(wave_power_low_high_TIP(:));
[ii,jj] = find(wave_power_low_high_TIP==mx);
    
pos = [pos; a1+ii,b1+jj];
    
    
Positions_low_high_TIP{iii} = pos;


end

