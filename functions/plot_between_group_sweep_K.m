clc
close all
clear

%% Varying strength of pleiotropy;

% Load colours
colours;

%% Within-group parameters

two_cell_experiment = 0;
cooperation_experiment = 1;

y_0 = 0;
y_by = 1;
y_end = 600;

rho = 1;
K = [100, 200, 500];
s = 0.95;
phi = 0;

mu = 0.0001;
nu = 0.01;

num_traits = 3;
num_types = 2^num_traits;

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

%% Parameter sweeps.

lambda = [10:5:50];
phi = [0.0:0.1:1];

%% Get data.

focalVariable = K;

data = zeros(length(focalVariable), length(lambda), length(phi), num_traits);

for k = 1:length(focalVariable)

    for j = 1:length(lambda)

        for i = 1:length(phi)

            % Save the processed data.
            filename_between_output = sprintf('group_results/lamg%g_Kg%d_sg%g_alphg%g_gamg%g_costg%g_two%d_coop%d_yend%d_rho%g_K%d_s%g_p%g_mu%g_nu%g', ...
                lambda(j), K_g, s_g, alpha, gamma, c_g, two_cell_experiment, cooperation_experiment, y_end, rho, K(k), s, phi(i), mu, nu);
            filename_between_output = strrep(filename_between_output, '.', '-');
            filename_between_output = strcat(filename_between_output, '.mat');

            if isfile([filename_between_output])

                file = load([filename_between_output], 'z_bar_t');
                data(k, j, i, :) = squeeze(file.z_bar_t(end, :));

            end

        end

    end

end

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
figname = sprintf('figure_between_group_K_sweep_two_%s_coop_%s_mu_%s_nu_%s_cost_%s_alpha_%s_gamma_%s', num2str(two_cell_experiment), num2str(cooperation_experiment), num2str(mu), num2str(nu), num2str(c_g), num2str(alpha), num2str(gamma));

% Set position % x_0 y_0 width heigth
fig = figure(1);
fig.Units = 'Centimeters';
fig.Position = [0, 0, width, height];

%% Subplot 1

ind = 1;
trait = 1;
valueK = 1;
z_bar_t = squeeze(data(valueK, :, :, trait));


% Convert linear index to row and column
[col, row] = ind2sub([numColPanels numRowPanels], ind);

% Positions.
posX = 1.5 + (col - 1) * panelSize;
posY = height - row * panelSize;

% Create the axis for the subplot [left bottom width height]
sp = subplot(numRowPanels, numColPanels, ind);
sp.Units = 'Centimeters';
sp.OuterPosition(1) = posX;
% sp.OuterPosition(2) = posY;
sp.OuterPosition(3) = panelSize;
sp.OuterPosition(4) = panelSize;
sp.Position(1) = sp.OuterPosition(1) + (panelSize - axisSize); % left plus some margin
sp.Position(2) = sp.OuterPosition(2) + (panelSize - axisSize); % bottom plus some margin
sp.Position(3) = axisSize;
sp.Position(4) = axisSize;

position = sp.Position;

h = heatmap(phi, flip(lambda), flip(z_bar_t, 1), 'GridVisible', 'off', 'CellLabelColor', 'none', 'FontSize', 8);
h.ColorbarVisible = 'off';
h.Colormap = squeeze(V(trait, :, :));
caxis(h, [0, 1])
h.XDisplayLabels = repmat({''}, size(h.XData));  %remove row labels

position = h.Position;

ax = axes();
ax.Units = 'Centimeters';
ax.Position = position;
ax.Visible = 'off';

plot(ax, [0.1 0.1], [0 1], ':', 'Color', hex2rgb('#212121'))
xlim([0 1])
set(gca, 'visible', 'off')

% fig_label = "A";
% text(0 - 0.45, 1.00, fig_label, 'Units', 'normalized', 'FontWeight', 'b', 'FontSize', 10)

title = "Cooperative";
text(0.5, 1.11, title, 'Units', 'normalized', 'FontSize', 8, 'HorizontalAlignment','center')

%% Subplot 2

ind = 2;
trait = 2;
valueK = 1;
z_bar_t = squeeze(data(valueK, :, :, trait));


% Convert linear index to row and column
[col, row] = ind2sub([numColPanels numRowPanels], ind);

% Positions.
posX = 1.5 + (col - 1) * panelSize;
posY = height - row * panelSize;

% Create the axis for the subplot [left bottom width height]
sp = subplot(numRowPanels, numColPanels, ind);
sp.Units = 'Centimeters';
sp.OuterPosition(1) = posX;
% sp.OuterPosition(2) = posY;
sp.OuterPosition(3) = panelSize;
sp.OuterPosition(4) = panelSize;
sp.Position(1) = sp.OuterPosition(1) + (panelSize - axisSize); % left plus some margin
sp.Position(2) = sp.OuterPosition(2) + (panelSize - axisSize); % bottom plus some margin
sp.Position(3) = axisSize;
sp.Position(4) = axisSize;

position = sp.Position;

h = heatmap(phi, flip(lambda), flip(z_bar_t, 1), 'GridVisible', 'off', 'CellLabelColor', 'none', 'FontSize', 8);
h.ColorbarVisible = 'off';
h.Colormap = squeeze(V(trait, :, :));
caxis(h, [0, 1])
h.XDisplayLabels = repmat({''}, size(h.XData));  %remove row labels
h.YDisplayLabels = repmat({''}, size(h.YData));  %remove y-labels

position = h.Position;

ax = axes();
ax.Units = 'Centimeters';
ax.Position = position;
ax.Visible = 'off';

plot(ax, [0.1 0.1], [0 1], ':', 'Color', hex2rgb('#212121'))
xlim([0 1])
set(gca, 'visible', 'off')

title = "Trait";
text(0.5, 1.31, title, 'Units', 'normalized', 'FontSize', 10, 'HorizontalAlignment','center','FontWeight', 'b')

title = "Private";
text(0.5, 1.11, title, 'Units', 'normalized', 'FontSize', 8, 'HorizontalAlignment','center')

%% Subplot 3

ind = 3;
trait = 3;
valueK = 1;
z_bar_t = squeeze(data(valueK, :, :, trait));


% Convert linear index to row and column
[col, row] = ind2sub([numColPanels numRowPanels], ind);

% Positions.
posX = 1.5 + (col - 1) * panelSize;
posY = height - row * panelSize;

% Create the axis for the subplot [left bottom width height]
sp = subplot(numRowPanels, numColPanels, ind);
sp.Units = 'Centimeters';
sp.OuterPosition(1) = posX;
% sp.OuterPosition(2) = posY;
sp.OuterPosition(3) = panelSize;
sp.OuterPosition(4) = panelSize;
sp.Position(1) = sp.OuterPosition(1) + (panelSize - axisSize); % left plus some margin
sp.Position(2) = sp.OuterPosition(2) + (panelSize - axisSize); % bottom plus some margin
sp.Position(3) = axisSize;
sp.Position(4) = axisSize;

position = sp.Position;

h = heatmap(phi, flip(lambda), flip(z_bar_t, 1), 'GridVisible', 'off', 'CellLabelColor', 'none', 'FontSize', 8);
h.ColorbarVisible = 'off';
h.Colormap = squeeze(V(trait, :, :));
caxis(h, [0, 1])
h.XDisplayLabels = repmat({''}, size(h.XData));  %remove row labels
h.YDisplayLabels = repmat({''}, size(h.YData));  %remove y-labels

position = h.Position;

ax = axes();
ax.Units = 'Centimeters';
ax.Position = position;
ax.Visible = 'off';

plot(ax, [0.1 0.1], [0 1], ':', 'Color', hex2rgb('#212121'))
xlim([0 1])
set(gca, 'visible', 'off')

title = "Pleiotropy";
text(0.5, 1.11, title, 'Units', 'normalized', 'FontSize', 8, 'HorizontalAlignment','center')

fig_label = "\it K = 100";
text(1 + 0.11, 0.5, fig_label, 'Units', 'normalized', 'FontSize', 8, 'Rotation', 270, 'HorizontalAlignment', 'center')

%% Subplot 4

ind = 4;
trait = 1;
valueK = 2;
z_bar_t = squeeze(data(valueK, :, :, trait));


% Convert linear index to row and column
[col, row] = ind2sub([numColPanels numRowPanels], ind);

% Positions.
posX = 1.5 + (col - 1) * panelSize;
posY = height - 1.32 * row * panelSize;

% Create the axis for the subplot [left bottom width height]
sp = subplot(numRowPanels, numColPanels, ind);
sp.Units = 'Centimeters';
sp.OuterPosition(1) = posX;
sp.OuterPosition(2) = posY;
sp.OuterPosition(3) = panelSize;
sp.OuterPosition(4) = panelSize;
sp.Position(1) = sp.OuterPosition(1) + (panelSize - axisSize); % left plus some margin
sp.Position(2) = sp.OuterPosition(2) + (panelSize - axisSize); % bottom plus some margin
sp.Position(3) = axisSize;
sp.Position(4) = axisSize;

position = sp.Position;

h = heatmap(phi, flip(lambda), flip(z_bar_t, 1), 'GridVisible', 'off', 'CellLabelColor', 'none', 'FontSize', 8);
h.ColorbarVisible = 'off';
h.Colormap = squeeze(V(trait, :, :));
caxis(h, [0, 1])
h.XDisplayLabels = repmat({''}, size(h.XData));  %remove row labels

position = h.Position;

ax = axes();
ax.Units = 'Centimeters';
ax.Position = position;
ax.Visible = 'off';

plot(ax, [0.1 0.1], [0 1], ':', 'Color', hex2rgb('#212121'))
xlim([0 1])
set(gca, 'visible', 'off')

text([-0.27], 0.11, 'Group lifespan, \lambda', 'Rotation', 90, 'VerticalAlignment', 'middle', 'FontSize', 8, 'FontName', 'Helvetica');


% fig_label = "B";
% text(0 - 0.45, 1.00, fig_label, 'Units', 'normalized', 'FontWeight', 'b', 'FontSize', 10)

%% Subplot 5

ind = 5;
trait = 2;
valueK = 2;
z_bar_t = squeeze(data(valueK, :, :, trait));

% Convert linear index to row and column
[col, row] = ind2sub([numColPanels numRowPanels], ind);

% Positions.
posX = 1.5 + (col - 1) * panelSize;
posY = height - 1.32 * row * panelSize;

% Create the axis for the subplot [left bottom width height]
sp = subplot(numRowPanels, numColPanels, ind);
sp.Units = 'Centimeters';
sp.OuterPosition(1) = posX;
sp.OuterPosition(2) = posY;
sp.OuterPosition(3) = panelSize;
sp.OuterPosition(4) = panelSize;
sp.Position(1) = sp.OuterPosition(1) + (panelSize - axisSize); % left plus some margin
sp.Position(2) = sp.OuterPosition(2) + (panelSize - axisSize); % bottom plus some margin
sp.Position(3) = axisSize;
sp.Position(4) = axisSize;

position = sp.Position;

h = heatmap(phi, flip(lambda), flip(z_bar_t, 1), 'GridVisible', 'off', 'CellLabelColor', 'none', 'FontSize', 8);
h.ColorbarVisible = 'off';
h.Colormap = squeeze(V(trait, :, :));
caxis(h, [0, 1])
h.XDisplayLabels = repmat({''}, size(h.XData));  %remove row labels
h.YDisplayLabels = repmat({''}, size(h.YData));  %remove y-labels

ax = axes();
ax.Units = 'Centimeters';
ax.Position = position;
ax.Visible = 'off';

plot(ax, [0.1 0.1], [0 1], ':', 'Color', hex2rgb('#212121'))
xlim([0 1])
set(gca, 'visible', 'off')

%% Subplot 6

ind = 6;
trait = 3;
valueK = 2;
z_bar_t = squeeze(data(valueK, :, :, trait));

% Convert linear index to row and column
[col, row] = ind2sub([numColPanels numRowPanels], ind);

% Positions.
posX = 1.5 + (col - 1) * panelSize;
posY = height - 1.32 * row * panelSize;

% Create the axis for the subplot [left bottom width height]
sp = subplot(numRowPanels, numColPanels, ind);
sp.Units = 'Centimeters';
sp.OuterPosition(1) = posX;
sp.OuterPosition(2) = posY;
sp.OuterPosition(3) = panelSize;
sp.OuterPosition(4) = panelSize;
sp.Position(1) = sp.OuterPosition(1) + (panelSize - axisSize); % left plus some margin
sp.Position(2) = sp.OuterPosition(2) + (panelSize - axisSize); % bottom plus some margin
sp.Position(3) = axisSize;
sp.Position(4) = axisSize;

position = sp.Position;

h = heatmap(phi, flip(lambda), flip(z_bar_t, 1), 'GridVisible', 'off', 'CellLabelColor', 'none', 'FontSize', 8);
h.ColorbarVisible = 'off';
h.Colormap = squeeze(V(trait, :, :));
caxis(h, [0, 1])
h.XDisplayLabels = repmat({''}, size(h.XData));  %remove row labels
h.YDisplayLabels = repmat({''}, size(h.YData));  %remove y-labels

ax = axes();
ax.Units = 'Centimeters';
ax.Position = position;
ax.Visible = 'off';

plot(ax, [0.1 0.1], [0 1], ':', 'Color', hex2rgb('#212121'))
xlim([0 1])
set(gca, 'visible', 'off')

fig_label = "\it K = 200";
text(1 + 0.11, 0.5, fig_label, 'Units', 'normalized', 'FontSize', 8, 'Rotation', 270, 'HorizontalAlignment', 'center')

fig_label = "Group size";
text(1 + 0.31, 0.5, fig_label, 'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'b', 'Rotation', 270, 'HorizontalAlignment', 'center')

%% Subplot 7

ind = 7;
trait = 1;
valueK = 3;
z_bar_t = squeeze(data(valueK, :, :, trait));

% Convert linear index to row and column
[col, row] = ind2sub([numColPanels numRowPanels], ind);

% Positions.
posX = 1.5 + (col - 1) * panelSize;
posY = height - 1.22 * row * panelSize;

% Create the axis for the subplot [left bottom width height]
sp = subplot(numRowPanels, numColPanels, ind);
sp.Units = 'Centimeters';
sp.OuterPosition(1) = posX;
sp.OuterPosition(2) = posY;
sp.OuterPosition(3) = panelSize;
sp.OuterPosition(4) = panelSize;
sp.Position(1) = sp.OuterPosition(1) + (panelSize - axisSize); % left plus some margin
sp.Position(2) = sp.OuterPosition(2) + (panelSize - axisSize); % bottom plus some margin
sp.Position(3) = axisSize;
sp.Position(4) = axisSize;

position = sp.Position;

h = heatmap(phi, flip(lambda), flip(z_bar_t, 1), 'GridVisible', 'off', 'CellLabelColor', 'none', 'FontSize', 8);
h.ColorbarVisible = 'off';
h.Colormap = squeeze(V(trait, :, :));
caxis(h, [0, 1])
% h.XDisplayLabels = repmat({''}, size(h.XData));  %remove x-labels
% h.YDisplayLabels = repmat({''}, size(h.YData));  %remove y-labels

position = h.Position;

ax = axes();
ax.Units = 'Centimeters';
ax.Position = position;
ax.Visible = 'off';

plot(ax, [0.1 0.1], [0 1], ':', 'Color', hex2rgb('#212121'))
xlim([0 1])
set(gca, 'visible', 'off')

% fig_label = "C";
% text(0 - 0.45, 1.00, fig_label, 'Units', 'normalized', 'FontWeight', 'b', 'FontSize', 10)

%% Colorbar 7

ax = axes();
ax.Units = 'Centimeters';
ax.Position = position;
ax.Visible = 'off';

c = colorbar(ax);
c.Ticks = 0:0.2:1;
c.Units = 'Centimeters';
c.Location = 'southoutside'; % Horizontal
c.Position(2) = c.Position(2) - 2.5;
colormap(ax, squeeze(V(trait, :, :)))

%% Subplot 8

ind = 8;
trait = 2;
valueK = 3;
z_bar_t = squeeze(data(valueK, :, :, trait));

% Convert linear index to row and column
[col, row] = ind2sub([numColPanels numRowPanels], ind);

% Positions.
posX = 1.5 + (col - 1) * panelSize;
posY = height - 1.22 * row * panelSize;

% Create the axis for the subplot [left bottom width height]
sp = subplot(numRowPanels, numColPanels, ind);
sp.Units = 'Centimeters';
sp.OuterPosition(1) = posX;
sp.OuterPosition(2) = posY;
sp.OuterPosition(3) = panelSize;
sp.OuterPosition(4) = panelSize;
sp.Position(1) = sp.OuterPosition(1) + (panelSize - axisSize); % left plus some margin
sp.Position(2) = sp.OuterPosition(2) + (panelSize - axisSize); % bottom plus some margin
sp.Position(3) = axisSize;
sp.Position(4) = axisSize;

position = sp.Position;

h = heatmap(phi, flip(lambda), flip(z_bar_t, 1), 'GridVisible', 'off', 'CellLabelColor', 'none', 'FontSize', 8);
h.ColorbarVisible = 'off';
h.Colormap = squeeze(V(trait, :, :));
caxis(h, [0, 1])
h.YDisplayLabels = repmat({''}, size(h.YData));  %remove row labels
% h.XDisplayLabels = repmat({''}, size(h.XData));  %remove row labels

xlabel('Pleiotropy strength, \phi')

ax = axes();
ax.Units = 'Centimeters';
ax.Position = position;
ax.Visible = 'off';

plot(ax, [0.1 0.1], [0 1], ':', 'Color', hex2rgb('#212121'))
xlim([0 1])
set(gca, 'visible', 'off')

%% Colorbar 8

ax = axes();
ax.Units = 'Centimeters';
ax.Position = position;
ax.Visible = 'off';

c = colorbar(ax);
c.Ticks = 0:0.2:1;
c.Units = 'Centimeters';
c.Location = 'southoutside'; % Horizontal
c.Position(2) = c.Position(2) - 2.5;
colormap(ax, squeeze(V(trait, :, :)))

%% Subplot 9

ind = 9;
trait = 3;
valueK = 3;
z_bar_t = squeeze(data(valueK, :, :, trait));

% Convert linear index to row and column
[col, row] = ind2sub([numColPanels numRowPanels], ind);

% Positions.
posX = 1.5 + (col - 1) * panelSize;
posY = height - 1.22 * row * panelSize;

% Create the axis for the subplot [left bottom width height]
sp = subplot(numRowPanels, numColPanels, ind);
sp.Units = 'Centimeters';
sp.OuterPosition(1) = posX;
sp.OuterPosition(2) = posY;
sp.OuterPosition(3) = panelSize;
sp.OuterPosition(4) = panelSize;
sp.Position(1) = sp.OuterPosition(1) + (panelSize - axisSize); % left plus some margin
sp.Position(2) = sp.OuterPosition(2) + (panelSize - axisSize); % bottom plus some margin
sp.Position(3) = axisSize;
sp.Position(4) = axisSize;

position = sp.Position;

h = heatmap(phi, flip(lambda), flip(z_bar_t, 1), 'GridVisible', 'off', 'CellLabelColor', 'none', 'FontSize', 8);
h.ColorbarVisible = 'off';
h.Colormap = squeeze(V(trait, :, :));
caxis(h, [0, 1])
% h.XDisplayLabels = repmat({''}, size(h.XData));  %remove row labels
h.YDisplayLabels = repmat({''}, size(h.YData));  %remove row labels

ax = axes();
ax.Units = 'Centimeters';
ax.Position = position;
ax.Visible = 'off';

plot(ax, [0.1 0.1], [0 1], ':', 'Color', hex2rgb('#212121'))
xlim([0 1])
set(gca, 'visible', 'off')

fig_label = "\it K = 500";
text(1 + 0.11, 0.5, fig_label, 'Units', 'normalized', 'FontSize', 8, 'Rotation', 270, 'HorizontalAlignment', 'center')

%% Colorbar 9

ax = axes();
ax.Units = 'Centimeters';
ax.Position = position;
ax.Visible = 'off';

c = colorbar(ax);
c.Ticks = 0:0.2:1;
c.Units = 'Centimeters';
c.Location = 'southoutside'; % Horizontal
c.Position(2) = c.Position(2) - 2.5;
colormap(ax, squeeze(V(trait, :, :)))

%% Export

export_fig(['figures/between-group/', figname], '-transparent', '-tiff', '-r300')
% export_fig(['figures/between-group/', figname], '-transparent', '-eps', '-r300')
