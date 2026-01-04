
clc

data_dam = load('./data/AL0625_1_dam.mat');
data_undam = load('./data/AL0625_1_undam.mat');


%%

wave_data_ = data_dam.AL0625_1_dam;
wave_data_dam = reshape(wave_data_,sqrt(size(wave_data_,1)),sqrt(size(wave_data_,1)),size(wave_data_,2));

wave_data_ = data_undam.AL0625_1_undam;
wave_data_undam = reshape(wave_data_,sqrt(size(wave_data_,1)),sqrt(size(wave_data_,1)),size(wave_data_,2));


wave_power_dam = sum(wave_data_dam.^2,3);
wave_power_undam = sum(wave_data_undam.^2,3);

%%

mx = max(wave_power_dam(:));
[ii,jj] = find(wave_power_dam==mx);

wave_data_dam = wave_data_dam(:,:,1:100);

%%

lx = linspace(0,99,size(wave_data_dam,1));
ly = linspace(0,99,size(wave_data_dam,2));

axisFontsize = 15;

%%

figure,
set(gcf, 'Units', 'centimeters', 'Position', [10, 10, 13, 10])
imagesc(lx,ly,wave_data_dam(:,:,5)/max(wave_data_dam(:,:,5),[],'all'));
xlabel('\fontsize{16}Plate width [mm]')
ylabel('\fontsize{16}Plate length [mm]')
axis xy;
axis square;
ax = gca;
%colorbar;
ax.FontSize=axisFontsize;
exportgraphics(gcf,'images/wave_at_5_dam.png','Resolution',600) 

figure,
set(gcf, 'Units', 'centimeters', 'Position', [10, 10, 13, 10])
imagesc(lx,ly,wave_data_dam(:,:,10)/max(wave_data_dam(:,:,10),[],'all'));
xlabel('\fontsize{16}Plate width [mm]')
ylabel('\fontsize{16}Plate length [mm]')
axis xy;
axis square;
ax = gca;
%colorbar;
ax.FontSize=axisFontsize;
exportgraphics(gcf,'images/wave_at_10_dam.png','Resolution',600) 

figure,
set(gcf, 'Units', 'centimeters', 'Position', [10, 10, 13, 10])
imagesc(lx,ly,wave_data_dam(:,:,13)/max(wave_data_dam(:,:,13),[],'all'));
xlabel('\fontsize{16}Plate width [mm]')
ylabel('\fontsize{16}Plate length [mm]')
axis xy;
axis square;
ax = gca;
%colorbar;
ax.FontSize=axisFontsize;
exportgraphics(gcf,'images/wave_at_13_dam.png','Resolution',600) 


figure,
set(gcf, 'Units', 'centimeters', 'Position', [10, 10, 13, 10])
imagesc(lx,ly,wave_power_dam/max(wave_power_dam(:)));
xlabel('\fontsize{16}Plate width [mm]')
ylabel('\fontsize{16}Plate length [mm]')
axis xy;
axis square;
ax = gca;
colorbar;
ax.FontSize=axisFontsize;
exportgraphics(gcf,'images/wave_energy_dam.png','Resolution',600) 




%%

figure,
set(gcf, 'Units', 'centimeters', 'Position', [10, 10, 13, 10])
imagesc(lx,ly,wave_data_undam(:,:,5)/max(wave_data_undam(:,:,5),[],'all'));
xlabel('\fontsize{16}Plate width [mm]')
ylabel('\fontsize{16}Plate length [mm]')
axis xy;
axis square;
ax = gca;
%colorbar;
ax.FontSize=axisFontsize;
exportgraphics(gcf,'images/wave_at_5_undam.png','Resolution',600) 

figure,
set(gcf, 'Units', 'centimeters', 'Position', [10, 10, 13, 10])
imagesc(lx,ly,wave_data_undam(:,:,10)/max(wave_data_undam(:,:,10),[],'all'));
xlabel('\fontsize{16}Plate width [mm]')
ylabel('\fontsize{16}Plate length [mm]')
axis xy;
axis square;
ax = gca;
%colorbar;
ax.FontSize=axisFontsize;
exportgraphics(gcf,'images/wave_at_10_undam.png','Resolution',600) 

figure,
set(gcf, 'Units', 'centimeters', 'Position', [10, 10, 13, 10])
imagesc(lx,ly,wave_data_undam(:,:,13)/max(wave_data_dam(:,:,13),[],'all'));
xlabel('\fontsize{16}Plate width [mm]')
ylabel('\fontsize{16}Plate length [mm]')
axis xy;
axis square;
ax = gca;
%colorbar;
ax.FontSize=axisFontsize;
exportgraphics(gcf,'images/wave_at_13_undam.png','Resolution',600) 

figure,
set(gcf, 'Units', 'centimeters', 'Position', [10, 10, 13, 10])
imagesc(lx,ly,wave_power_undam/max(wave_power_undam(:)));
xlabel('\fontsize{16}Plate width [mm]')
ylabel('\fontsize{16}Plate length [mm]')
axis xy;
axis square;
ax = gca;
colorbar;
ax.FontSize=axisFontsize;
exportgraphics(gcf,'images/wave_energy_undam.png','Resolution',600) 



%%
