
clc

addpath(genpath( './all_libs/WIMF'))

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

Lx = 4*del2(eye(size(wave_data_dam,2)));
Ly = 4*del2(eye(size(wave_data_dam,1)));
Lt = 4*del2(eye(size(wave_data_dam,3)));

[V_x,Dx] = eig(Lx);
[V_y,Dy] = eig(Ly);
[V_t,Dt] = eig(Lt);

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



%% Primary Source Localization Images

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

% Define subplot labels
labels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'};

% Define the 100x100 square boundaries
square_size = 100;

% Create a background matrix A 
A = wave_data_dam(:,:,5); 

% Load all results files except secondary_results.mat
result_files = dir('./results/*_results.mat');
result_files = result_files(~strcmp({result_files.name}, 'secondary_results.mat'));

% Sort files for consistent ordering
[~, sort_idx] = sort({result_files.name});
result_files = result_files(sort_idx);

% Define colors for rectangles and background (matching reference image)
rect_color = [0.53, 0.66, 0.87];    % Soft blue for inside rectangle
bg_color = [0.82, 0.9, 0.68];       % Pale greenish yellow for outside region
rect_alpha = 0.85;                   % Higher opacity for rectangle
bg_alpha = 0.85;                     % Higher opacity for background region
matrix_alpha = 0.10;                 % Lower transparency for matrix A overlay

% Define different marker styles
marker_styles = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*'};

% Create output directory if it doesn't exist
output_dir = './saved_plots';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Loop through each of the 6 rectangles
for indx = 1:6
    % Create a new figure for each plot
    fig = figure('Position', [100, 100, 360, 360], 'Visible', 'off'); 
    hold on;
    
    % Get the 4 vertices of the current rectangle
    vertices = squeeze(centers(:, :, indx)); % 4x2 matrix
    
    % Extract coordinates
    x_coords = vertices(:, 2);
    y_coords = vertices(:, 1);
    rect_x = [x_coords; x_coords(1)]; % Close the rectangle
    rect_y = [y_coords; y_coords(1)];
    
    % FIRST: Fill the background region (outside rectangle)
    h_bg = fill([0, square_size, square_size, 0], [0, 0, square_size, square_size], ...
         bg_color, 'FaceAlpha', bg_alpha, 'EdgeColor', 'none');
    set(h_bg, 'HandleVisibility', 'off'); % Exclude from legend
    
    % SECOND: Fill the rectangle on top
    h_rect = fill(rect_x, rect_y, rect_color, 'FaceAlpha', rect_alpha, ...
         'EdgeColor', 'none');
    set(h_rect, 'HandleVisibility', 'off'); % Exclude from legend
    
    % THIRD: Display background matrix A with transparency ON TOP
    h_img = imagesc([0 square_size], [0 square_size], A);
    set(h_img, 'AlphaData', matrix_alpha);  % Make matrix A transparent
    colormap(gca, 'parula');  % Using parula colormap
    set(h_img, 'HandleVisibility', 'off'); % Exclude from legend
    
    % Load and plot positions from results files
    for file_idx = 1:length(result_files)
        % Load the file
        data = load(fullfile('./results', result_files(file_idx).name));
        
        % Extract the Positions cell array
        if isfield(data, 'Positions')
            positions_cell = data.Positions;
            
            % Get the positions for this index
            if indx <= length(positions_cell)
                positions = positions_cell{indx}; % n x 2 matrix
                
                if ~isempty(positions)
                    % Extract filename
                    filename = result_files(file_idx).name;
                    
                    % % Special handling for primary_secondary_results.mat
                    % if strcmp(filename, 'primary_secondary_results.mat')
                    %     label_base = 'WIRL';
                    % else
                        % Extract filename without extension and drop everything after last _
                        %filename_no_ext = filename(1:end-12); % Remove '_results.mat'
                        underscore_pos = strfind(filename, '_');
                        if ~isempty(underscore_pos)
                            label_base = filename(1:underscore_pos(end)-1);
                        else
                            label_base = filename;
                        end
                    % end
                    
                    % Get marker style for this file
                    marker_style = marker_styles{mod(file_idx-1, length(marker_styles)) + 1};
                    
                    % Plot each row of positions
                    n_positions = size(positions, 1);
                    
                    if n_positions > 1
                        % Use same marker with different colors for multiple positions
                        row_colors = lines(n_positions);
                        for row_idx = 1:n_positions
                            if row_idx == n_positions
                               label_name = sprintf('%s-full', label_base);
                            else
                                label_name = sprintf('%s_%d', label_base, row_idx);
                            end
                            
                            % Plot the position with same marker but different color (no edge)
                            plot(positions(row_idx, 2), positions(row_idx, 1), ...
                                 marker_style, 'MarkerSize', 8, ...
                                 'Color', row_colors(row_idx, :), ...
                                 'MarkerFaceColor', row_colors(row_idx, :), ...
                                 'MarkerEdgeColor', 'none', ...
                                 'DisplayName', label_name);
                        end
                    else
                        % Single position - use file-specific color
                        position_colors = lines(length(result_files));
                        label_name = label_base;
                        label_name = strrep(label_name, '_', '-');

                        % Plot the position (no edge)
                        plot(positions(1, 2), positions(1, 1), ...
                             marker_style, 'MarkerSize', 8, ...
                             'Color', position_colors(file_idx, :), ...
                             'MarkerFaceColor', position_colors(file_idx, :), ...
                             'MarkerEdgeColor', 'none', ...
                             'DisplayName', label_name);
                    end
                end
            end
        end
    end
    
    % Set axis properties
    axis square;
    axis xy;
    xlim([0, square_size]);
    ylim([0, square_size]);
    xlabel('Plate width [mm]', 'FontSize', 22);
    ylabel('Plate length [mm]', 'FontSize', 22);
    title(sprintf('Mask %s', labels{indx}), 'FontSize', 16);
    
    % text(50, -15, sprintf('Mask %s', labels{indx}), ...
    %  'FontSize', 16, ...
    %  'HorizontalAlignment', 'center');

    grid on;
    ax = gca;  % Get current axes
    ax.GridAlpha = 0.9;  % Set grid opacity (0 = transparent, 1 = opaque)
    ax.GridColor = [0, 0, 0];  % Black grid
    %ax.GridLineStyle = '--';     % Solid line (default is '-')
    ax.FontSize = 13;
    % Only show legend on the first plot
    if indx == 1
        leg = legend('Location', 'southwest', 'FontSize', 8);
        leg.Box = 'on';
        leg.Color = [1 1 1 0.8]; % White background with slight transparency
    end
    
    hold off;
    
    % Save the figure
    output_filename = fullfile(output_dir, sprintf('mask_%s.png', labels{indx}));
    saveas(fig, output_filename);
    
    % Optional: Save as high-resolution image
    print(fig, fullfile(output_dir, sprintf('mask_%s_highres.png', labels{indx})), ...
          '-dpng', '-r300');
    
    % Close the figure
    close(fig);
    
    fprintf('Saved plot %d to %s\n', indx, output_filename);
end

fprintf('All plots saved successfully to %s\n', output_dir);


%%

% Load all result files
all_files = dir('./results/*_results.mat');



% Files to exclude
exclude_files = {'WIRL_results.mat'};

% Filter out excluded files
valid_files = {};
for i = 1:length(all_files)
    if ~ismember(all_files(i).name, exclude_files)
        valid_files{end+1} = all_files(i).name;
    end
end

% Specify index (1 to 6)
indx = 3;  % Change this as needed
slice_idx = 10;

% Output directory
output_dir = './recon_figures';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

for i = 1:length(valid_files)
    filename = valid_files{i};
    filepath = fullfile('./results', filename);
    data = load(filepath);
    
    if strcmp(filename, 'WIMF_recon_results.mat')
        recon = inv_spectral_wave_transform(data.all_Ds{indx} * data.all_xs{indx}, V_x, V_y, V_t);
    else
        recon = data.all_recons{indx};
    end
    
    % Create figure (invisible)
    fig = figure('Visible', 'off', 'Position', [100, 100, 400, 350]);
    
    imagesc(recon(:,:,slice_idx));
    axis square;
    axis xy;
    colorbar;
    colormap('parula');
    xlabel('Plate Width [mm]');
    ylabel('Plate Length [mm]');
    
    % Clean up title
    [~, name, ~] = fileparts(filename);
    name = strrep(name, '_results', '');
    name = strrep(name, 'WIMF_recon', 'WIRL');
    title(strrep(name, '_', '-'), 'FontSize', 11);
    
    % Save as PNG (high resolution for IEEE)
    output_filename = fullfile(output_dir, sprintf('%s_indx%d_slice%d.png', name, indx, slice_idx));
    exportgraphics(fig, output_filename, 'Resolution', 300);
    
    close(fig);
    fprintf('Saved: %s\n', output_filename);
end

fprintf('Done! All figures saved to %s\n', output_dir);


%%

% Load all defect files
all_files = dir('./results/*_defect.mat');


% Mask Calc

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

% Slice index for 3D arrays
slice = 5;

unmask = ones(size(wave_data_dam,1),size(wave_data_dam,2));

a1 = centers(1,1,slice);
a4 = centers(4,1,slice);
b1 = centers(1,2,slice);
b2 = centers(2,2,slice);


unmask(a1:a4,b1:b2,:) = 0;

% Overlay color (e.g., white)
overlay_color = [1, 1, 1];  % RGB: white
overlay_alpha = 0.5;  % transparency (0 = invisible, 1 = opaque)

% Output directory
output_dir = './defect_figures';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

indx = 3;

for i = 1:length(all_files)
    filename = all_files(i).name;
    filepath = fullfile('./results', filename);
    data = load(filepath);
    
    % Get the data
    disp_data = data.WavePower_est{slice};
    
    % Check dimensions
    if ndims(disp_data) == 2
        plot_data = disp_data;
    else
        plot_data = disp_data(:,:,slice);
    end
    
    % Create figure (invisible)
    fig = figure('Visible', 'off', 'Position', [100, 100, 400, 350]);
    
    imagesc(plot_data);
    axis xy;
    axis square;
    hold on;
    
    % Create RGB overlay
    overlay = cat(3, overlay_color(1)*ones(size(unmask)), ...
                     overlay_color(2)*ones(size(unmask)), ...
                     overlay_color(3)*ones(size(unmask)));
    h = image(overlay);
    set(h, 'AlphaData', overlay_alpha * unmask);
    
    hold off;
    colorbar;
    colormap('jet');
    xlabel('Plate Width [mm]');
    ylabel('Plate Length [mm]');
    
    % Clean up title
    [~, name, ~] = fileparts(filename);
    name = strrep(name, '_defect', '');
    title(strrep(name, '_', '-'), 'FontSize', 11);
    
    % Save as PNG
    output_filename = fullfile(output_dir, sprintf('%s_defect_slice%d.png', name, slice));
    exportgraphics(fig, output_filename, 'Resolution', 300);
    
    close(fig);
    fprintf('Saved: %s\n', output_filename);
end

fprintf('Done! All figures saved to %s\n', output_dir);

%%

rmpath(genpath( './all_libs/WIMF'))