clc
clear

data_dam = load('./data/AL0625_1_dam.mat');
data_undam = load('./data/AL0625_1_undam.mat');

rng(42);

%% Preprocess data and mask to have the required shape for further steps

% Wave data with damage present
wave_data_dam_ = data_dam.AL0625_1_dam;
wave_data_dam = reshape(wave_data_dam_,sqrt(size(wave_data_dam_,1)),sqrt(size(wave_data_dam_,1)),size(wave_data_dam_,2));

wave_data_dam = wave_data_dam(:,:,1:100);

% Wave data without damage

wave_data_undam_ = data_undam.AL0625_1_undam;
wave_data_undam = reshape(wave_data_undam_,sqrt(size(wave_data_undam_,1)),sqrt(size(wave_data_undam_,1)),size(wave_data_undam_,2));

wave_data_undam = wave_data_undam(:,:,1:100);


% Estimate Primary source and damage location using energy method
wave_power = sum(wave_data_dam.^2,3);

mx = max(wave_power(:));
[trueI,trueJ] = find(wave_power==mx);

mg = min(wave_power(:));
[trueI_,trueJ_] = find(wave_power==mg);

Lx = 4*del2(eye(size(wave_data_dam,2)));
Ly = 4*del2(eye(size(wave_data_dam,1)));
Lt = 4*del2(eye(size(wave_data_dam,3)));

[V_x,Dx] = eig(Lx);
[V_y,Dy] = eig(Ly);
[V_t,Dt] = eig(Lt);

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

%%
indx = 5;

a1 = centers(1,1,indx); b1 = centers(1,2,indx);
a2 = centers(2,1,indx); b2 = centers(2,2,indx);
a3 = centers(3,1,indx); b3 = centers(3,2,indx);
a4 = centers(4,1,indx); b4 = centers(4,2,indx);

lx = linspace(0,99,size(wave_data_undam,1));
ly = linspace(0,99,size(wave_data_undam,2));

axisFontsize = 15;

Trans = ones(size(wave_data_undam));
Trans(a1:a4,b1:b2,:) = 0;

set(gcf, 'Units', 'centimeters', 'Position', [10, 5, 13*2, 10*2])
hold on;
imagesc(lx,ly,wave_data_undam(:,:,5)/max(wave_data_undam(:,:,5),[],'all'),'AlphaData',0.3);
imagesc(lx,ly,Trans(:,:,5)/max(Trans(:,:,5),[],'all'),'AlphaData',0.3);
xlim([0,99])
ylim([0,99])
xlabel('\fontsize{16}Plate width [mm]')
ylabel('\fontsize{16}Plate length [mm]')
axis xy;
axis square;
ax = gca;
%colorbar;
grid on;
ax.FontSize=axisFontsize;

%%
addpath(genpath( './all_libs/WIMF'))

wave_data_dam_ = wave_data_dam;
wave_data_undam_ = wave_data_undam;

wave_data_dam_(a1:a4,b1:b2,:) = 0;
wave_data_undam_(a1:a4,b1:b2,:) = 0;

Trans = ones(size(wave_data_dam));
Trans(a1:a4,b1:b2,:) = 0;


[Dc_dam,x_dam,C_dam] = waveInformedMatFac(wave_data_dam_,Trans,V_x,V_y,V_t,'count',3,'gradient_descent',true);
[Dc_undam,x_undam,C_undam] = waveInformedMatFac(wave_data_undam_,Trans,V_x,V_y,V_t,'count',3,'gradient_descent',true);

pos_defect = [];
[C_undam_permuted, perm_idx] = permute_diagonal_nearest(C_dam, C_undam);

Dc_undam = Dc_undam(:,perm_idx);
x_undam = x_undam(perm_idx);

all_D_dam{indx} = Dc_dam;
all_D_undam{indx} = Dc_undam;

all_x_dam{indx} = x_dam;
all_x_undam{indx} = x_undam;

all_C_dam{indx} = C_dam;
all_C_undam{indx} = C_undam;

for jjj=1:size(Dc_dam,2)

    power_D_dam = inv_spectral_wave_transform(Dc_dam(:,jjj)*x_dam(jjj),V_x,V_y,V_t);
    power_D_undam = inv_spectral_wave_transform(Dc_undam(:,jjj)*x_undam(jjj),V_x,V_y,V_t);

    power_D = power_D_dam - power_D_undam;

    wave_power = sum(power_D(a1:a4,b1:b2,8:13).^2,3);
    mx = max(wave_power(:));
    [ii,jj] = find(wave_power==mx);

    pos_defect = [pos_defect; a1+ii,b1+jj];

    % mg = min(wave_power(:));
    % [ij,ji] = find(wave_power==mg);
    % pos_defect = [pos_defect; a1+ij,b1+ji];

end

reconstructedData_dam = inv_spectral_wave_transform(Dc_dam*x_dam,V_x,V_y,V_t);
reconstructedData_undam = inv_spectral_wave_transform(Dc_undam*x_undam,V_x,V_y,V_t);

defect = reconstructedData_dam - reconstructedData_undam;

wave_power = sum(defect(a1:a4,b1:b2,8:12).^2,3);
% wave_power = abs(sum(reconstructedData_dam(a1:a4,b1:b2,1:15).^2,3)-sum(reconstructedData_undam(a1:a4,b1:b2,1:15).^2,3));
mx = max(wave_power(:));
[ii,jj] = find(wave_power==mx);
pos_defect = [pos_defect; a1+ii,b1+jj];
    
Positions_defect{indx} = pos_defect;

rmpath(genpath( './all_libs/WIMF'))

