function [] = subplot_within_dynamics_n(ind, figWidth, figHeight, panelSize, axisSize, numColPanels, numRowPanels, n_y, idx_25, num_types, y_lim, x_lim, y_lab, x_lab, title, fig_label, y_tick, y_tick_label, x_tick, x_tick_label)

% Load colours
colours;

% Convert linear index to row and column
[col, row] = ind2sub([numColPanels numRowPanels], ind);

% Positions.
posX = 1.5 + (col - 1) * panelSize;
posY = figHeight - 0.5 - row * panelSize;

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

% Create the plot
for i = 1:num_types
    pl = plot(n_y(:, i), '-', 'Color', genotype_colour{i}, 'LineWidth', 1.1);
    hold on
end

hold on;
if fig_label ~= ""
    if idx_25 > 3
        xline(idx_25, '-.', '25% Invasion', 'FontSize', 8);
    end
    
else
    if idx_25 > 3
        xline(idx_25, '-.', 'FontSize', 8);
    end
    
end
hold off;

ylim(y_lim)
xlim(x_lim)

if y_lab ~= ""
    lh = ylabel(y_lab, 'FontSize', 10);
    lh.Position(1) = lh.Position(1) - 0.05;
end

if x_lab ~= ""
    xlabel(x_lab, 'FontSize', 10)
end

if fig_label ~= ""
    text(0 - 0.55, 1.00, fig_label, 'Units', 'normalized', 'FontWeight', 'b', 'FontSize', 10)
end

if title ~= ""
    text(0.35, 1.1, title, 'Units', 'normalized', 'FontSize', 8)
end

% Modify the axes
ax = gca;
ax.FontSize = 8; % Sets the tick label font size
ax.YTick = y_tick;
ax.YTickLabel = y_tick_label;
ax.XTick = x_tick;
ax.XTickLabel = x_tick_label;

end

