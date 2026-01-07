clc
clear

data = load('./data/AL0625_1_dam.mat');


%% Preprocess data and mask to have the required shape for further steps

wave_data_ = data.AL0625_1_dam;
wave_data = reshape(wave_data_,sqrt(size(wave_data_,1)),sqrt(size(wave_data_,1)),size(wave_data_,2));


wave_power = sum(wave_data.^2,3);

wave_data = wave_data(:,:,1:100);

mx = max(wave_power(:));
[trueI,trueJ] = find(wave_power==mx);

mg = min(wave_power(:));
[trueI_,trueJ_] = find(wave_power==mg);

Lx = 4*del2(eye(size(wave_data,2)));
Ly = 4*del2(eye(size(wave_data,1)));
Lt = 4*del2(eye(size(wave_data,3)));

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

for iii = 1:6

    addpath(genpath( './all_libs/WIMF'))

    a1 = centers(1,1,iii); b1 = centers(1,2,iii);
    a2 = centers(2,1,iii); b2 = centers(2,2,iii);
    a3 = centers(3,1,iii); b3 = centers(3,2,iii);
    a4 = centers(4,1,iii); b4 = centers(4,2,iii);

    wave_data_ = wave_data;
    wave_data_(a1:a4,b1:b2,:) = 0;
    
    Trans = ones(size(wave_data));
    Trans(a1:a4,b1:b2,:) = 0;
    
    [Dc,x,C] = waveInformedMatFac(wave_data_,Trans,V_x,V_y,V_t,'count',4,'gradient_descent',true,'tolerance',0.00001);
    
    all_Ds{iii} = Dc;
    all_xs{iii} = x;
    all_Cs{iii} = C;
    %all_losses{iii} = l_h;

    pos_primary = [];
    pos_secondary = [];
    for jjj=1:size(Dc,2)
        
        power_D = inv_spectral_wave_transform(Dc(:,jjj),V_x,V_y,V_t);
        wave_power = sum(power_D(a1:a4,b1:b2,1:5).^2,3);
        mx = max(wave_power(:));
        [ii,jj] = find(wave_power==mx);
        pos_primary = [pos_primary; a1+ii,b1+jj];

        mg = min(wave_power(:));
        [ij,ji] = find(wave_power==mg);

        pos_secondary = [pos_secondary; a1+ij,b1+ji];
    end
    
    reconstructedData = inv_spectral_wave_transform(Dc*x,V_x,V_y,V_t);
    wave_power = sum(reconstructedData(a1:a4,b1:b2,1:5).^2,3);
    mx = max(wave_power(:));
    [ii,jj] = find(wave_power==mx);

    mg = min(wave_power(:));
    [ij,ji] = find(wave_power==mg);
    Positions_secondary{iii} = pos_secondary;
    
    all_recons{iii} = reconstructedData;

    pos_primary = [pos_primary; a1+ii,b1+jj];
    pos_secondary = [pos_secondary; a1+ij,b1+ji];
    
    Positions{iii} = pos_primary;
    
    
    psnr_vals{iii} = calculate_psnr(wave_data,reconstructedData);
    
    rmpath(genpath( './all_libs/WIMF'))
end

save('./results/primary_secondary_results.mat','Positions', ...
    'Positions_secondary',"psnr_vals",'all_Ds','all_xs','all_Cs');
