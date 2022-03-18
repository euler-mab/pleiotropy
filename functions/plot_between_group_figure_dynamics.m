clc
close all
clear

%% Varying strength of pleiotropy;

% Load colours
colours;

two_cell_experiment = 0;
cooperation_experiment = 1;

y_0 = 0;
y_by = 1;
y_end = 600;

rho = 1;
K = 500;
s = 0.95;
phi = [0.0, 0.5, 1.0];

mu = 0.0001;
nu = 0.01;

num_traits = 3; % number of binary traits
num_types = 2^num_traits; % number of cell types resulting from binary traits

Z = [
    0, 0, 0;
    0, 0, 1;
    0, 1, 0;
    0, 1, 1;
    1, 0, 0;
    1, 0, 1;
    1, 1, 0;
    1, 1, 1
    ];

% between-group

num_group_traits = 5;
num_group_types = num_types;
if two_cell_experiment == 1
    num_group_types = num_types * (num_types - 1) / 2 + num_types;
end

lambda = [10:5:50]; % parameter lambda in the text, affecting expected lifespan of groups
K_g = 1; % carrying capacity of group popularion (normalised to 1 in this workflow)
s_g = 0.95; % strength of between-group selection
c_g = 0.0; % cost of pleiotropy
gamma = 0; % parameter gamma in the text, affecting germ line vs. somatic cell
alpha = 0.0; % parameter alpha in the text, affecting age of reproductive maturity.

%%

for i = 1:length(lambda)

    % Store the outputs.
    data_n_group_t_y_output = {};
    data_z_bar_group_t_output = {};
    data_n_t_output = {};
    data_z_bar_t_output = {};

    for j = 1:length(phi)

        % load the between-group dynamics results
        filename_between_output = sprintf('group_results/lamg%g_Kg%d_sg%g_alphg%g_gamg%g_costg%g_two%d_coop%d_yend%d_rho%g_K%d_s%g_p%g_mu%g_nu%g', ...
            lambda(i), K_g, s_g, alpha, gamma, c_g, two_cell_experiment, cooperation_experiment, y_end, rho, K, s, phi(j), mu, nu);
        filename_between_output = strrep(filename_between_output, '.', '-');
        filename_between_output = strcat(filename_between_output, '.mat');

        if isfile(filename_between_output)

            file = load(filename_between_output, 'n_group_t_y', 'z_bar_group_t', 'n_t', 'z_bar_t');
            data_n_group_t_y_output{j} = file.n_group_t_y;
            data_z_bar_group_t_output{j} = file.z_bar_group_t;
            data_n_t_output{j} = file.n_t;
            data_z_bar_t_output{j} = file.z_bar_t;

        else

            % Do nothing

        end

    end


    %% Figure parameters

    x_end = 10000;

    %% Plot

    % Set width and height from Plos Biology.
    maxPanels = 6;
    panelSize = 19.05 / (1.0 * maxPanels);
    axisSize = 2.75;

    numRowPanels = 3;
    numColPanels = length(phi);

    width = 2 * panelSize * numColPanels;
    height = 1.5 * panelSize * numRowPanels;

    % Figure name
    figname = sprintf('figure_between_group_dynamics_lamg_%s_K_%s_coop_%s_c_%s', num2str(lambda(i)), num2str(K), num2str(cooperation_experiment), num2str(c_g));

    % Set position % x_0 y_0 width heigth
    fig = figure(1);
    fig.Units = 'Centimeters';
    fig.Position = [0, 0, width, height];

    %% Subplot 1

    ind = 1;
    x_t = data_n_t_output{1} ./ sum(data_n_t_output{1}, 2);

    x_end = 10000;

    y_lim = [0 1.1 * 1];
    x_lim = [0 x_end];
    y_lab = 'Frequency,\it x^c(t)';
    x_lab = '';
    title = '\phi = 0.0';
    fig_label = 'A';
    y_tick = 0:0.2:1;
    y_tick_label = 0:0.2:1;
    x_tick = [];
    x_tick_label = [];
    stretch = 0;
    subplot_between_dynamics_x(ind, width, height, panelSize, stretch, axisSize, numColPanels, numRowPanels, x_t, num_types, y_lim, x_lim, y_lab, x_lab, title, fig_label, y_tick, y_tick_label, x_tick, x_tick_label)

    %% Subplot 2

    ind = 2;
    x_t = data_n_t_output{2} ./ sum(data_n_t_output{2}, 2);

    x_end = 5000;

    y_lim = [0 1.1 * 1];
    x_lim = [0 x_end];
    y_lab = '';
    x_lab = '';
    title = '\phi = 0.5';
    fig_label = '';
    y_tick = [];
    y_tick_label = [];
    x_tick = [];
    x_tick_label = [];
    stretch = 1;
    subplot_between_dynamics_x(ind, width, height, panelSize, stretch, axisSize, numColPanels, numRowPanels, x_t, num_types, y_lim, x_lim, y_lab, x_lab, title, fig_label, y_tick, y_tick_label, x_tick, x_tick_label)

    %% Subplot 3

    ind = 3;
    x_t = data_n_t_output{3} ./ sum(data_n_t_output{3}, 2);

    x_end = 5000;

    y_lim = [0 1.1 * 1];
    x_lim = [0 x_end];
    y_lab = '';
    x_lab = '';
    title = '\phi = 1.0';
    fig_label = '';
    y_tick = [];
    y_tick_label = [];
    x_tick = [];
    x_tick_label = [];
    stretch = 1;
    subplot_between_dynamics_x(ind, width, height, panelSize, stretch, axisSize, numColPanels, numRowPanels, x_t, num_types, y_lim, x_lim, y_lab, x_lab, title, fig_label, y_tick, y_tick_label, x_tick, x_tick_label)



    %% Subplot 4

    ind = 4;
    z_bar_t = data_z_bar_t_output{1};
    x_t = data_n_t_output{1} ./ sum(data_n_t_output{1}, 2);

    x_end = 10000;

    y_lim = [0 1.1 * 1];
    x_lim = [0 x_end];
    y_lab = ['Avg trait,\it z^c(t)'];
    x_lab = '';
    title = '';
    fig_label = 'B';
    y_tick = 0:0.2:1;
    y_tick_label = 0:0.2:1;
    x_tick = [];
    x_tick_label = [];
    stretch = 0;
    subplot_between_dynamics_z_bar(ind, width, height, panelSize, stretch, axisSize, numColPanels, numRowPanels, z_bar_t, num_traits, y_lim, x_lim, y_lab, x_lab, title, fig_label, y_tick, y_tick_label, x_tick, x_tick_label)

    %% Subplot 5

    ind = 5;
    z_bar_t = data_z_bar_t_output{2};
    x_t = data_n_t_output{2} ./ sum(data_n_t_output{2}, 2);

    x_end = 5000;

    y_lim = [0 1.1 * 1];
    x_lim = [0 x_end];
    y_lab = '';
    x_lab = '';
    title = '';
    fig_label = '';
    y_tick = [];
    y_tick_label = [];
    x_tick = [];
    x_tick_label = [];
    stretch = 1;
    subplot_between_dynamics_z_bar(ind, width, height, panelSize, stretch, axisSize, numColPanels, numRowPanels, z_bar_t, num_traits, y_lim, x_lim, y_lab, x_lab, title, fig_label, y_tick, y_tick_label, x_tick, x_tick_label)

    %% Subplot 6

    ind = 6;
    z_bar_t = data_z_bar_t_output{3};
    x_t = data_n_t_output{3} ./ sum(data_n_t_output{3}, 2);

    x_end = 5000;

    y_lim = [0 1.1 * 1];
    x_lim = [0 x_end];
    y_lab = '';
    x_lab = '';
    title = '';
    fig_label = '';
    y_tick = [];
    y_tick_label = [];
    x_tick = [];
    x_tick_label = [];
    stretch = 1;
    subplot_between_dynamics_z_bar(ind, width, height, panelSize, stretch, axisSize, numColPanels, numRowPanels, z_bar_t, num_traits, y_lim, x_lim, y_lab, x_lab, title, fig_label, y_tick, y_tick_label, x_tick, x_tick_label)



    %% Subplot 7

    ind = 7;

    z_bar_group_t = data_z_bar_group_t_output{1};
    z_bar_t = data_z_bar_t_output{1};

    if size(z_bar_t) ~= [0 0]
        z_bar_mutation_load = z_bar_group_t(:, 1:3) - z_bar_t;
    else
        z_bar_mutation_load = 0 * z_bar_t;
    end

    x_t = data_n_t_output{1} ./ sum(data_n_t_output{1}, 2);

    x_end = 10000;

    y_lim = [-0.4 * 0.1 0.1 * 0.4];
    x_lim = [0 x_end];
    y_lab = 'Lifetime trait change';
    x_lab = 'Time,\it t';
    title = '';
    fig_label = 'C';
    y_tick = [-0.4 * 0.1 : 0.02 : 0.1 * 0.4];
    y_tick_label =  [-0.4 * 0.1 : 0.02 : 0.1 * 0.4];
    x_tick = 0:x_end/4:x_end;
    x_tick_label = 0:x_end/4:x_end;
    stretch = 0;
    subplot_between_dynamics_z_bar_load(ind, width, height, panelSize, stretch, axisSize, numColPanels, numRowPanels, z_bar_mutation_load, num_traits, y_lim, x_lim, y_lab, x_lab, title, fig_label, y_tick, y_tick_label, x_tick, x_tick_label)

    %% Subplot 8

    ind = 8;
    z_bar_group_t = data_z_bar_group_t_output{2};
    z_bar_t = data_z_bar_t_output{2};

    if size(z_bar_t) ~= [0 0]
        z_bar_mutation_load = z_bar_group_t(:, 1:3) - z_bar_t;
    else
        z_bar_mutation_load = 0 * z_bar_t;
    end

    x_t = data_n_t_output{2} ./ sum(data_n_t_output{2}, 2);

    x_end = 5000;

    y_lim = [-0.4 * 0.1 0.1 * 0.4];
    x_lim = [0 x_end];
    y_lab = '';
    x_lab = 'Time,\it t';
    title = '';
    fig_label = '';
    y_tick = [];
    y_tick_label =  [];
    x_tick = 0:x_end/4:x_end;
    x_tick_label = 0:x_end/4:x_end;
    stretch = 1;
    subplot_between_dynamics_z_bar_load(ind, width, height, panelSize, stretch, axisSize, numColPanels, numRowPanels, z_bar_mutation_load, num_traits, y_lim, x_lim, y_lab, x_lab, title, fig_label, y_tick, y_tick_label, x_tick, x_tick_label)

    %% Subplot 9

    ind = 9;

    z_bar_group_t = data_z_bar_group_t_output{3};
    z_bar_t = data_z_bar_t_output{3};

    if size(z_bar_t) ~= [0 0]
        z_bar_mutation_load = z_bar_group_t(:, 1:3) - z_bar_t;
    else
        z_bar_mutation_load = 0 * z_bar_t;
    end

    x_t = data_n_t_output{3} ./ sum(data_n_t_output{3}, 2);

    x_end = 5000;

    y_lim = [-0.4 * 0.1 0.1 * 0.4];
    x_lim = [0 x_end];
    y_lab = '';
    x_lab = 'Time,\it t';
    title = '';
    fig_label = '';
    y_tick = [];
    y_tick_label =  [];
    x_tick = 0:x_end/4:x_end;
    x_tick_label = 0:x_end/4:x_end;
    stretch = 1;
    subplot_between_dynamics_z_bar_load(ind, width, height, panelSize, stretch, axisSize, numColPanels, numRowPanels, z_bar_mutation_load, num_traits, y_lim, x_lim, y_lab, x_lab, title, fig_label, y_tick, y_tick_label, x_tick, x_tick_label)


    %% Export

    export_fig(['figures/between-group/', figname], '-transparent', '-tiff', '-r300')
    % export_fig(['figures/between-group/', figname], '-transparent', '-eps', '-r300')

end