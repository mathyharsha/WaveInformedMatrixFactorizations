clear all
close all
clc

%% ====================== Load data ==============================
clc

data = load('/Users/harshat/Downloads/UF_work/all_libs/AL0625_1_dam.mat');


%%

wave_data_ = data.AL0625_1_dam;
wave_data = reshape(wave_data_,sqrt(size(wave_data_,1)),sqrt(size(wave_data_,1)),size(wave_data_,2));

wave_data = wave_data(:,:,1:150);

wave_data_un = wave_data;

Trans = ones(size(wave_data));
Trans(6:90,20:90,:) = 0;

wave_data(6:90,20:90,:) = 0;

trans = reshape(Trans,size(Trans,1)*size(Trans,2)*size(Trans,3),1);

T = wave_data;

%%


X               =       T/max(abs(T(:)))                   ;
R               =       zeros(size(X))                     ;  % for multi tube noise
for i = 1:(size(X,3)/10)
    temp        =       (rand(size(X,1),size(X,2))>0.99)   ;
    R(:,:,i*10-9:i*10) = repmat(temp,[1 1 10])             ;
end
Z               =       randn(size(X))                     ;
Z               =       Z.*R                               ;
Xn              =       X+Z                                ;

lambda          =       1/sqrt(size(X,1))                  ;

%% ================ main process of rPCA =======================
[ L , S ]       =       tensor_rpca( X , lambda)          ;

%% ===================== Plot Result ===========================

fh = figure;
for i = 1:size(X,3) 
    subplot(221)
    imagesc(Xn(:,:,i));title('Noisy Video');
    colormap(gray);
    
    subplot(222)
    imagesc(1-Z(:,:,i));title('True Sparse Noise');
    colormap(gray);
    
    subplot(223)
    imagesc(L(:,:,i));title('Low Rank Component');
    colormap(gray);
      
    subplot(224)
    imagesc(1-S(:,:,i));title('Sparse Component')
    colormap(gray);
    
    pause(0.1);

end