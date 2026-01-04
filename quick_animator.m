
data_dam = load('./data/AL0625_1_dam.mat');
data_undam = load('./data/AL0625_1_undam.mat');

%% Load Data if Needed for some visualization

% Wave data with damage present
wave_data_dam_ = data_dam.AL0625_1_dam;
wave_data_dam = reshape(wave_data_dam_,sqrt(size(wave_data_dam_,1)),sqrt(size(wave_data_dam_,1)),size(wave_data_dam_,2));

wave_data_dam = wave_data_dam(:,:,1:100);

% Wave data without damage

wave_data_undam_ = data_undam.AL0625_1_undam;
wave_data_undam = reshape(wave_data_undam_,sqrt(size(wave_data_undam_,1)),sqrt(size(wave_data_undam_,1)),size(wave_data_undam_,2));

wave_data_undam = wave_data_undam(:,:,1:100);

%% Animate actual data
close all
figure('Position', [100, 100, 1200, 500]) 
for i = 1:floor(size(wave_data_undam,3))
    % Actual With damage
    subplot(1,3,1)
    imagesc(wave_data_dam(:,:,i))
    axis equal tight
    colorbar
    title(['Actual data - With Damage ' num2str(i)])
    
    % Actual Without damage
    subplot(1,3,2)
    imagesc(wave_data_undam(:,:,i))
    axis equal tight
    colorbar
    title(['Actual data - Without Damage ' num2str(i)])
    
    subplot(1,3,3)
    imagesc(wave_data_dam(:,:,i)-wave_data_undam(:,:,i))
    axis equal tight
    colorbar
    title(['The difference ' num2str(i)])

    drawnow  % Force update of figure
    pause(0.5)
end


%%




close all
figure('Position', [100, 100, 1200, 500])  % Create figure with specific size [left, bottom, width, height]

for i = 1:25 %floor(size(defect,3))
    % WIMF Reconstruction with damage
    imagesc(reconstructedData(a1:a4,b1:b2,i))
    axis equal tight
    colorbar
    title(['WIMF Reconstruction - With Damage ' num2str(i)])
    drawnow
    pause(0.1)
end




%% Animate WIMF Reconstruction

addpath(genpath( './all_libs/WIMF'))
chs_indx = 5;
ind = 2;
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

a1 = centers(1,1,chs_indx); b1 = centers(1,2,chs_indx);
a2 = centers(2,1,chs_indx); b2 = centers(2,2,chs_indx);
a3 = centers(3,1,chs_indx); b3 = centers(3,2,chs_indx);
a4 = centers(4,1,chs_indx); b4 = centers(4,2,chs_indx);

reconstructedData_dam = inv_spectral_wave_transform(all_D_dam{chs_indx}*all_x_dam{chs_indx},V_x,V_y,V_t); 
reconstructedData_undam = inv_spectral_wave_transform(all_D_undam{chs_indx}*all_x_undam{chs_indx},V_x,V_y,V_t);

close all
figure('Position', [100, 100, 1200, 500])  % Create figure with specific size [left, bottom, width, height]

for i = 1:25 %floor(size(defect,3))
    % WIMF Reconstruction with damage
    subplot(1,3,1)
    imagesc(reconstructedData_dam(a1:a4,b1:b2,i))
    axis equal tight
    colorbar
    title(['WIMF Reconstruction - With Damage ' num2str(i)])
    
    % Without damage
    subplot(1,3,2)
    imagesc(reconstructedData_undam(a1:a4,b1:b2,i))
    axis equal tight
    colorbar
    title(['WIMF Reconstruction - Without Damage ' num2str(i)])
    
    % Baseline subtracted
    subplot(1,3,3)
    imagesc(defect(a1:a4,b1:b2,i))
    axis equal tight
    colorbar
    title(['Baseline subtracted ' num2str(i)])
    drawnow  % Force update of figure
    pause(0.7)
end

rmpath(genpath( './all_libs/WIMF'))

% %% Load Data if Needed for some visualization
% 
% % Wave data with damage present
% wave_data_dam_ = data_dam.AL0625_1_dam;
% wave_data_dam = reshape(wave_data_dam_,sqrt(size(wave_data_dam_,1)),sqrt(size(wave_data_dam_,1)),size(wave_data_dam_,2));
% 
% wave_data_dam = wave_data_dam(:,:,1:100);
% 
% % Wave data without damage
% 
% wave_data_undam_ = data_undam.AL0625_1_undam;
% wave_data_undam = reshape(wave_data_undam_,sqrt(size(wave_data_undam_,1)),sqrt(size(wave_data_undam_,1)),size(wave_data_undam_,2));
% 
% wave_data_undam = wave_data_undam(:,:,1:100);

