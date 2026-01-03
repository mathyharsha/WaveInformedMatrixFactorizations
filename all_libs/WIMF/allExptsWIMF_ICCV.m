clc
clear

data = load('AL0625_1_dam.mat');


%% Preprocess data and mask to have the required shape for further steps

wave_data_ = data.AL0625_1_dam;
wave_data = reshape(wave_data_,sqrt(size(wave_data_,1)),sqrt(size(wave_data_,1)),size(wave_data_,2));

wave_data = wave_data(:,:,1:100);

wave_power = sum(wave_data.^2,3);

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

centers = [1;1;1;1]*[trueI,trueJ];
%%

for iii = 1:3

    [a1,b1] = mids(1,1,centers(1,1),centers(1,2));
    [a2,b2] = mids(1,100,centers(2,1),centers(2,2));
    [a3,b3] = mids(100,100,centers(3,1),centers(3,2));
    [a4,b4] = mids(100,1,centers(4,1),centers(4,2));
    
    centers(1,1) = a1; centers(1,2) = b1;
    centers(2,1) = a2; centers(2,2) = b2;
    centers(3,1) = a3; centers(3,2) = b3;
    centers(4,1) = a4; centers(4,2) = b4;

    if iii<3
        continue
    end

    wave_data_ = wave_data;
    wave_data_(a1:a4,b1:b2,:) = 0;
    
    Trans = ones(size(wave_data));
    Trans(a1:a4,b1:b2,:) = 0;
    
    
    [Dc,x,C] = waveInformedMatFac(wave_data_,Trans,V_x,V_y,V_t,'count',5);
    
    all_Ds{iii} = Dc;
    pos_primary = [];
    pos_secondary = [];
    for jjj=1:size(Dc,2)
        
        power_D = inv_spectral_wave_transform(Dc(:,jjj),V_x,V_y,V_t);
        wave_power = sum(power_D(a1:a4,b1:b2,1:10).^2,3);
        mx = max(wave_power(:));
        [ii,jj] = find(wave_power==mx);
        pos_primary = [pos_primary; a1+ii,b1+jj];

        mg = min(wave_power(:));
        [ij,ji] = find(wave_power==mg);

        pos_secondary = [pos_secondary; a1+ij,b1+ji];

    end
    
    reconstructedData = inv_spectral_wave_transform(Dc*x,V_x,V_y,V_t);
    wave_power = sum(reconstructedData(a1:a4,b1:b2,5:10).^2,3);
    mx = max(wave_power(:));
    [ii,jj] = find(wave_power==mx);

    mg = min(wave_power(:));
    [ij,ji] = find(wave_power==mg);
        
    pos_primary = [pos_primary; a1+ii,b1+jj];
    pos_secondary = [pos_secondary; a1+ij,b1+ji];
        
    Positions_primary{iii} = pos_primary;
    Positions_secondary{iii} = pos_secondary;

end

