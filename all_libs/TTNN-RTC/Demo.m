%% Robust tensor completion using transformed tensor singular value decomposition
% by Guangjing Song  Michael K. Ng  Xiongjun Zhang, Numerical Linear Algebra with Applications,
%  27(3):e2299, 2020.
% sGS-ADMM for robust tensor completion
% If you have any question, please contact to Xiongjun Zhang
% (xjzhang@mail.ccnu.edu.cn)
%%

clear all

addpath('tensor_toolbox_2.6');
randn('seed',2013); randn('seed',2013);

%%  Video data
% true data 
load('CarphoneG')
XT = CarphoneG./max(CarphoneG(:));

[n1 n2 n3] = size(XT);

SL = [0.1];  %% sparse ratio \gamma

for jj = 1:length(SL)

    
STS = sptenrand([n1 n2 n3], SL(jj));
STS = tensor(STS);

ZZ = XT + STS;
ZZ = ZZ.data;

% % XT = XT/n3;

dim = [n1, n2,n3];

%% initial parameters
beta = 0.05;


gamma = 1.618;
MaxIte = 500;
tol = 5e-4;


%%  Set  Omega
sr = 0.6;  %% ssampling ratio \rho
fprintf('Sampling ratio = %0.8e\n',sr);
fprintf('sparse corruption ratio = %0.8e\n',SL(jj));
temp = randperm(n1*n2*n3);
kks = round((sr)*n1*n2*n3);
mark = zeros(n1,n2,n3); 
mark(temp(1:kks)) = 1;

mu = 2/sqrt(max(n1,n2)*n3);


%% Initial point

X0 = zeros(n1,n2,n3);
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
    [X E k eta aaa] = TNN(ZZ,mark,dim,opts);
   toc

%% print the relative error, psnr
    Error = norm(X(:)-XT(:))/norm(XT(:));
    fprintf('Relative error = %0.8e\n',Error);
    PSNR = psnr(XT(:),X(:));
    fprintf('PSNR = %0.8e\n',PSNR);
    

%% U Transform

fprintf('===== t-SVD by unitary transform =====\n');

%% U matrix
O = tenmat(X,[3]); % unfolding
O = O.data;
[U D V] = svd(O,'econ');

mu = 25/sqrt(max(n1,n2)*n3);
opts.mu = mu;
%%
  tic
    [XU E k] = UTNN(U, ZZ,mark,dim,opts);
  toc

%% print the relative error, psnr

    Error = norm(XU(:)-XT(:))/norm(XT(:));
    fprintf('Relative error = %0.8e\n',Error);
    PSNR = psnr(XT(:),XU(:));
    fprintf('PSNR = %0.8e\n',PSNR);
 
end
 


