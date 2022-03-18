
K = 100;

% Store the outputs.
data_n_y_output = {};
data_N_y_output = {};
data_x_y_output = {};
data_z_bar_y_output = {};

data_fe_b_output = {};
data_fe_count_output = {};
data_fe_d_output = {};

for i = 1:length(phi)

    % Load the abundances - verbose.
    [ n_y, N_y, x_y, z_bar_y, fe_b_output, fe_d_output, fe_count_output  ] = load_processed_within_group_dynamics(two_cell_experiment, cooperation_experiment, t_end, rho, K, s_c, phi(i), mu, nu);

    data_n_y_output{i} = n_y;
    data_N_y_output{i} = N_y;
    data_x_y_output{i} = x_y;
    data_z_bar_y_output{i} = z_bar_y;
    
    data_fe_b_output{i} = fe_b_output;    
    data_fe_d_output{i} = fe_d_output;
    data_fe_count_output{i} = fe_count_output;

end

%%

for k = 1:num_types

%% Figure parameters

focal_type = k;

%% Plot

% Set width and height from Plos Biology.
maxPanels = 6;
panelSize = 19.05 / (1.02 * maxPanels);
axisSize = 2.75;

numRowPanels = 4;
numColPanels = 3;

width = 1.5 * panelSize * numColPanels;
height = 1.5 * panelSize * numRowPanels;

% Figure name
figname = sprintf('figure_within_group_dynamics_type_%s_K_%s_coop_%s', num2str(focal_type), num2str(K), num2str(cooperation_experiment));

% Set position % x_0 y_0 width heigth
fig = figure(1);
fig.Units = 'Centimeters';
fig.Position(3) = width; % width
fig.Position(4) = height; % height
fig.PaperPositionMode = 'auto';




%% Subplot 1

ind = 1;
n_y = data_n_y_output{1};
n_y = squeeze(n_y(:, focal_type, :));
x_y = data_x_y_output{1};
x_y = squeeze(x_y(:, focal_type, :));

idx_25 = find(x_y(:, focal_type) < 0.75);
if length(idx_25) > 0
    idx_25 = idx_25(1);
else
    idx_25 = squeeze(size(x_y, 1));
end

y_lim = [0 1.1 * K];
x_lim = [1 1.1 * 50]; 
y_lab = ['Abundance,\it n^c(y)'];
x_lab = '';
title = '\phi = 0.0';
fig_label = 'A';
y_tick = 1:(K/4):1.1 * K;
y_tick_label = 0:(K/4):1.1 * K;
x_tick = [];
x_tick_label = [];
subplot_within_dynamics_n(ind, width, height, panelSize, axisSize, numColPanels, numRowPanels, n_y, idx_25, num_types, y_lim, x_lim, y_lab, x_lab, title, fig_label, y_tick, y_tick_label, x_tick, x_tick_label)


%% Subplot 2

ind = 2;
n_y = data_n_y_output{2};
n_y = squeeze(n_y(:, focal_type, :));
x_y = data_x_y_output{2};
x_y = squeeze(x_y(:, focal_type, :));

idx_25 = find(x_y(:, focal_type) < 0.75);
if length(idx_25) > 0
    idx_25 = idx_25(1);
else
    idx_25 = squeeze(size(x_y, 1));
end

y_lim = [0 1.1 * K];
x_lim = [1 1.1 * 50]; 
y_lab = '';
x_lab = '';
title = '\phi = 0.5';
fig_label = '';
y_tick = [];
y_tick_label = [];
x_tick = [];
x_tick_label = [];
subplot_within_dynamics_n(ind, width, height, panelSize, axisSize, numColPanels, numRowPanels, n_y, idx_25, num_types, y_lim, x_lim, y_lab, x_lab, title, fig_label, y_tick, y_tick_label, x_tick, x_tick_label)

%% Subplot 3

ind = 3;
n_y = data_n_y_output{3};
n_y = squeeze(n_y(:, focal_type, :));
x_y = data_x_y_output{3};
x_y = squeeze(x_y(:, focal_type, :));

idx_25 = find(x_y(:, focal_type) < 0.75);
if length(idx_25) > 0
    idx_25 = idx_25(1);
else
    idx_25 = squeeze(size(x_y, 1));
end

y_lim = [0 1.1 * K];
x_lim = [1 1.1 * 50]; 
y_lab = '';
x_lab = '';
title = '\phi = 1.0';
fig_label = '';
y_tick = [];
y_tick_label = [];
x_tick = [];
x_tick_label = [];
subplot_within_dynamics_n(ind, width, height, panelSize, axisSize, numColPanels, numRowPanels, n_y, idx_25, num_types, y_lim, x_lim, y_lab, x_lab, title, fig_label, y_tick, y_tick_label, x_tick, x_tick_label)





%% Subplot 4

ind = 4;
x_y = data_x_y_output{1};
x_y = squeeze(x_y(:, focal_type, :));

idx_25 = find(x_y(:, focal_type) < 0.75);
if length(idx_25) > 0
    idx_25 = idx_25(1);
else
    idx_25 = squeeze(size(x_y, 1));
end

y_lim = [0 1.1 * 1];
x_lim = [1 1.1 * 50]; 
y_lab = ['Frequency,\it x^c(y)'];
x_lab = '';
title = '';
fig_label = 'B';
y_tick = 0:0.2:1;
y_tick_label = 0:0.2:1;
x_tick = [];
x_tick_label = [];
subplot_within_dynamics_x(ind, width, height, panelSize, axisSize, numColPanels, numRowPanels, x_y, idx_25, num_types, y_lim, x_lim, y_lab, x_lab, title, fig_label, y_tick, y_tick_label, x_tick, x_tick_label)


%% Subplot 5

ind = 5;
x_y = data_x_y_output{2};
x_y = squeeze(x_y(:, focal_type, :));

idx_25 = find(x_y(:, focal_type) < 0.75);
if length(idx_25) > 0
    idx_25 = idx_25(1);
else
    idx_25 = squeeze(size(x_y, 1));
end

y_lim = [0 1.1 * 1];
x_lim = [1 1.1 * 50]; 
y_lab = '';
x_lab = '';
title = '';
fig_label = '';
y_tick = [];
y_tick_label = [];
x_tick = [];
x_tick_label = [];
subplot_within_dynamics_x(ind, width, height, panelSize, axisSize, numColPanels, numRowPanels, x_y, idx_25, num_types, y_lim, x_lim, y_lab, x_lab, title, fig_label, y_tick, y_tick_label, x_tick, x_tick_label)

%% Subplot 6

ind = 6;
x_y = data_x_y_output{3};
x_y = squeeze(x_y(:, focal_type, :));

idx_25 = find(x_y(:, focal_type) < 0.75);
if length(idx_25) > 0
    idx_25 = idx_25(1);
else
    idx_25 = squeeze(size(x_y, 1));
end

y_lim = [0 1.1 * 1];
x_lim = [1 1.1 * 50]; 
y_lab = '';
x_lab = '';
title = '';
fig_label = '';
y_tick = [];
y_tick_label = [];
x_tick = [];
x_tick_label = [];
subplot_within_dynamics_x(ind, width, height, panelSize, axisSize, numColPanels, numRowPanels, x_y, idx_25, num_types, y_lim, x_lim, y_lab, x_lab, title, fig_label, y_tick, y_tick_label, x_tick, x_tick_label)





%% Subplot 7

ind = 7;
z_bar_y = data_z_bar_y_output{1};
z_bar_y = squeeze(z_bar_y(:, focal_type, :));

x_y = data_x_y_output{1};
x_y = squeeze(x_y(:, focal_type, :));

idx_25 = find(x_y(:, focal_type) < 0.75);
if length(idx_25) > 0
    idx_25 = idx_25(1);
else
    idx_25 = squeeze(size(x_y, 1));
end

y_lim = [0 1.1 * 1];
x_lim = [1 1.1 * 50]; 
y_lab = ['Avg trait,\it z^c(y)'];
x_lab = '';
title = '';
fig_label = 'C';
y_tick = 0:0.2:1;
y_tick_label = 0:0.2:1;
x_tick = [1, 10:10:50];
x_tick_label = [0, 10:10:50];
subplot_within_dynamics_z_bar(ind, width, height, panelSize, axisSize, numColPanels, numRowPanels, z_bar_y, idx_25, num_traits, y_lim, x_lim, y_lab, x_lab, title, fig_label, y_tick, y_tick_label, x_tick, x_tick_label)


%% Subplot 8

ind = 8;
z_bar_y = data_z_bar_y_output{2};
z_bar_y = squeeze(z_bar_y(:, focal_type, :));

x_y = data_x_y_output{2};
x_y = squeeze(x_y(:, focal_type, :));

idx_25 = find(x_y(:, focal_type) < 0.75);
if length(idx_25) > 0
    idx_25 = idx_25(1);
else
    idx_25 = squeeze(size(x_y, 1));
end

y_lim = [0 1.1 * 1];
x_lim = [1 1.1 * 50]; 
y_lab = '';
x_lab = 'Group age,\it y';
title = '';
fig_label = '';
y_tick = [];
y_tick_label = [];
x_tick = [1, 10:10:50];
x_tick_label = [0, 10:10:50];
subplot_within_dynamics_z_bar(ind, width, height, panelSize, axisSize, numColPanels, numRowPanels, z_bar_y, idx_25, num_traits, y_lim, x_lim, y_lab, x_lab, title, fig_label, y_tick, y_tick_label, x_tick, x_tick_label)

%% Subplot 9

ind = 9;
z_bar_y = data_z_bar_y_output{3};
z_bar_y = squeeze(z_bar_y(:, focal_type, :));

x_y = data_x_y_output{3};
x_y = squeeze(x_y(:, focal_type, :));

idx_25 = find(x_y(:, focal_type) < 0.75);
if length(idx_25) > 0
    idx_25 = idx_25(1);
else
    idx_25 = squeeze(size(x_y, 1));
end

y_lim = [0 1.1 * 1];
x_lim = [1 1.1 * 50]; 
y_lab = '';
x_lab = '';
title = '';
fig_label = '';
y_tick = [];
y_tick_label = [];
x_tick = [1, 10:10:50];
x_tick_label = [0, 10:10:50];
subplot_within_dynamics_z_bar(ind, width, height, panelSize, axisSize, numColPanels, numRowPanels, z_bar_y, idx_25, num_traits, y_lim, x_lim, y_lab, x_lab, title, fig_label, y_tick, y_tick_label, x_tick, x_tick_label)

    
if k ~= 1

    %% Subplot 10

    ind = 10;

    fe_b = data_fe_b_output{1};    
    fe_d = data_fe_d_output{1};
    fe_count = data_fe_count_output{1};

    y_lim = [0 1.1 * 1];
    x_lim = [-6 6];
    y_lab = 'Mutation frequency';
    x_lab = '';
    title = '';
    fig_label = 'D';
    y_tick = 0:0.2:1;
    y_tick_label = 0:0.2:1;
    x_tick = [-6:3:-3, 0, 3:3:6];
    x_tick_label = [-6:3:-3, 0, 3:3:6];
    subplot_within_dynamics_dfe2(ind, width, height, panelSize, axisSize, numColPanels, numRowPanels, fe_b, fe_d, fe_count, Z, focal_type, num_traits, y_lim, x_lim, y_lab, x_lab, title, fig_label, y_tick, y_tick_label, x_tick, x_tick_label)

    %% Subplot 11

    ind = 11;

    fe_b = data_fe_b_output{2};    
    fe_d = data_fe_d_output{2};
    fe_count = data_fe_count_output{2};

    y_lim = [0 1.1 * 1];
    x_lim = [-6 6];
    y_lab = '';
    x_lab = 'Fitness effect (DFE)';
    title = '';
    fig_label = '';
    y_tick = [];
    y_tick_label = [];
    x_tick = [-6:3:-3, 0, 3:3:6];
    x_tick_label = [-6:3:-3, 0, 3:3:6];
    subplot_within_dynamics_dfe2(ind, width, height, panelSize, axisSize, numColPanels, numRowPanels, fe_b, fe_d, fe_count, Z, focal_type, num_traits, y_lim, x_lim, y_lab, x_lab, title, fig_label, y_tick, y_tick_label, x_tick, x_tick_label)

    %% Subplot 12

    ind = 12;

    fe_b = data_fe_b_output{3};    
    fe_d = data_fe_d_output{3};
    fe_count = data_fe_count_output{3};

    y_lim = [0 1.1 * 1];
    x_lim = [-6 6];
    y_lab = '';
    x_lab = '';
    title = '';
    fig_label = '';
    y_tick = [];
    y_tick_label = [];
    x_tick = [-6:3:-3, 0, 3:3:6];
    x_tick_label = [-6:3:-3, 0, 3:3:6];
    subplot_within_dynamics_dfe2(ind, width, height, panelSize, axisSize, numColPanels, numRowPanels, fe_b, fe_d, fe_count, Z, focal_type, num_traits, y_lim, x_lim, y_lab, x_lab, title, fig_label, y_tick, y_tick_label, x_tick, x_tick_label)

end
%% Export

export_fig(['figures/within-group/', figname], '-transparent', '-tiff', '-r300')

end