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

%% Computation of Eigenvalues and Eigenvectors of Second Derivative Matrix
%  to use for reconstruction

Lx = 4*del2(eye(size(wave_data_dam,2)));
Ly = 4*del2(eye(size(wave_data_dam,1)));
Lt = 4*del2(eye(size(wave_data_dam,3)));

[V_x,Dx] = eig(Lx);
[V_y,Dy] = eig(Ly);
[V_t,Dt] = eig(Lt);




%% Print summary
fprintf('\n========================================\n');
fprintf('Data Processing Complete\n');
fprintf('========================================\n');
fprintf('Wave data shape: %d x %d x %d\n', size(wave_data_dam,1), size(wave_data_dam,2), size(wave_data_dam,3));
fprintf('Primary source location: (%d, %d)\n', trueI, trueJ);
fprintf('Secondary source location: (%d, %d)\n', trueI_, trueJ_);
fprintf('========================================\n\n');



%%

% Load the .mat files
load('./results/primary_secondary_results.mat', 'Positions_primary');
load('./results/secondary_results.mat', 'Positions_defect');


% Define mask labels
mask_labels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'};

% Initialize cell arrays to store table data
nCases = 6; % Number of cases

% Preallocate cell arrays for the tables
primaryTableData = cell(nCases, 1);
defectTableData = cell(nCases, 1);

% Process Positions_primary data
for i = 1:nCases
    positions = Positions_primary{i};
    if ~isempty(positions)
        % Format as coordinate pairs: (x1, y1), (x2, y2), ...
        coordStrings = cell(size(positions, 1), 1);
        for j = 1:size(positions, 1)
            coordStrings{j} = sprintf('(%.2f, %.2f)', positions(j, 1), positions(j, 2));
        end
        primaryTableData{i} = strjoin(coordStrings, '; ');
    else
        primaryTableData{i} = 'No detections';
    end
end

% Process Positions_defect data
for i = 1:nCases
    positions = Positions_defect{i};
    if ~isempty(positions)
        % Format as coordinate pairs
        coordStrings = cell(size(positions, 1), 1);
        for j = 1:size(positions, 1)
            coordStrings{j} = sprintf('(%.2f, %.2f)', positions(j, 1), positions(j, 2));
        end
        defectTableData{i} = strjoin(coordStrings, '; ');
    else
        defectTableData{i} = 'No detections';
    end
end

% Create the tables using mask labels
% Table 1: Primary damage positions
primaryTable = table(mask_labels', primaryTableData, ...
    'VariableNames', {'Mask', 'Estimated_Primary_Source_Locations'});

% Table 2: Secondary defect positions
defectTable = table(mask_labels', defectTableData, ...
    'VariableNames', {'Mask', 'Locations_Secondary_Source_Locations'});

% Display the tables
disp('Table 1: Primary Damage Detection Results for WIMF (N_D = 4)');
disp(primaryTable);

fprintf('\n\n');

disp('Table 2: Secondary Defect Detection Results for WIMF (N_D = 1)');
disp(defectTable);

% Optional: Save tables to file
writetable(primaryTable, './results/primary_damage_table.csv');
writetable(defectTable, './results/secondary_defect_table.csv');

fprintf('\nTables saved to CSV files.\n');