clear all
close all
clc

%% ====================== Load data ==============================
clc

data = load('AL0625_1_dam.mat');


%%

wave_data_ = data.AL0625_1_dam;
wave_data = reshape(wave_data_,sqrt(size(wave_data_,1)),sqrt(size(wave_data_,1)),size(wave_data_,2));

wave_data = wave_data(:,:,1:15);

wave_data_un = wave_data;

Trans = ones(size(wave_data));
Trans(6:90,20:90,:) = 0;

wave_data(6:90,20:90,:) = 0;

trans = reshape(Trans,size(Trans,1)*size(Trans,2)*size(Trans,3),1);

normalize = max(wave_data(:));

b = wave_data(:)/max(wave_data(:));

%%

[n1,n2,n3]             =        size(wave_data)               ;
alpha                  =        1                             ;
maxItr                 =        1000                          ; % maximum iteration
rho                    =        0.01                          ;
myNorm                 =        'tSVD_1'                      ; % dont change for now

%% ================ main process of completion =======================
X   =    tensor_cpl_admm( trans , b , rho , alpha , ...
                     [n1,n2,n3] , maxItr , myNorm , 0 );
X                      =        X * normalize                 ;
X                      =        reshape(X,[n1,n2,n3])         ;
            
X_dif                  =        wave_data-X                           ;
RSE                    =        norm(X_dif(:))/norm( wave_data(:))     ;


%%

figure,
for i = 1:15
    subplot(221);imagesc(wave_data_un(6:90,20:90,i));axis off;
    title('Original Video');
    subplot(222);imagesc(wave_data(6:90,20:90,i)) ;axis off;
    title('Sampled Video');
    subplot(224);imagesc(X(6:90,20:90,i));axis off;
    title('Recovered Video');
    pause(.5);
end
            