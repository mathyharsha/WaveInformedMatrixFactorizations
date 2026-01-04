clc
clear

data_dam = load('./data/AL0625_1_dam.mat');
data_undam = load('./data/AL0625_1_undam.mat');

rng(42);

%% Preprocess data and mask to have the required shape for further steps

% Wave data with damage present
wave_data_dam_ = data_dam.AL0625_1_dam;
wave_data_dam = reshape(wave_data_dam_,sqrt(size(wave_data_dam_,1)),sqrt(size(wave_data_dam_,1)),size(wave_data_dam_,2));

% wave_data_dam = wave_data_dam(:,:,1:100);

% Wave data without damage

wave_data_undam_ = data_undam.AL0625_1_undam;
wave_data_undam = reshape(wave_data_undam_,sqrt(size(wave_data_undam_,1)),sqrt(size(wave_data_undam_,1)),size(wave_data_undam_,2));

% wave_data_undam = wave_data_undam(:,:,1:100);


% Estimate Primary source and damage location using energy method
wave_power = sum(wave_data_dam.^2,3);

mx = max(wave_power(:));
[trueI,trueJ] = find(wave_power==mx);

mg = min(wave_power(:));
[trueI_,trueJ_] = find(wave_power==mg);

%% Print summary
fprintf('\n========================================\n');
fprintf('Data Processing Complete\n');
fprintf('========================================\n');
fprintf('Wave data shape: %d x %d x %d\n', size(wave_data_dam,1), size(wave_data_dam,2), size(wave_data_dam,3));
fprintf('Primary source location: (%d, %d)\n', trueI, trueJ);
fprintf('Secondary source location: (%d, %d)\n', trueI_, trueJ_);
fprintf('========================================\n\n');


%% Print estimated defect positions for each mask

load('results/secondary_results.mat')
fprintf('========================================\n');
fprintf('Estimated Secondary Source Locations\n');
fprintf('========================================\n');
mask_labels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'};
for i = 1:length(Positions_defect)
    pos = Positions_defect{i};
    fprintf('Mask %s: (%d, %d)\n', mask_labels{i}, pos(1), pos(2));
end
fprintf('========================================\n\n');