function [] = subplot_within_dynamics_dfe(ind, figWidth, figHeight, panelSize, axisSize, numColPanels, numRowPanels, fe_b, fe_d, fe_count, Z, focal_type, num_traits, y_lim, x_lim, y_lab, x_lab, title, fig_label, y_tick, y_tick_label, x_tick, x_tick_label)

% Load colours
colours;

% Convert linear index to row and column
[col, row] = ind2sub([numColPanels numRowPanels], ind);

% Positions.
posX = 1.5 + (col - 1) * panelSize;
posY = figHeight - 1 - row * panelSize;

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




expected_age = 50;

fe_b_output_1 = fe_b;
fe_count_output_1 = fe_count;
fe_d_output_1 = fe_d;

fe_count_1 = sum(fe_count_output_1(1:expected_age, focal_type, :, Z(:, 1) == 1, ~Z(:, 1) == 1), 3);
fe_b_1 = sum(fe_b_output_1(1:expected_age, focal_type, :, Z(:, 1) == 1, ~Z(:, 1) == 1), 3);
fe_d_1 = sum(fe_d_output_1(1:expected_age, focal_type, :, Z(:, 1) == 1, ~Z(:, 1) == 1), 3);

fe_b_1 = fe_b_1(fe_count_1 ~= 0) ./ fe_count_1(fe_count_1 ~= 0);
fe_d_1 = fe_d_1(fe_count_1 ~= 0) ./ fe_count_1(fe_count_1 ~= 0);
fe_count_1 = fe_count_1(fe_count_1 ~= 0);

fe_count_2 = sum(fe_count_output_1(1:expected_age, focal_type, :, Z(:, 2) == 1, ~Z(:, 2) == 1), 3);
fe_b_2 = sum(fe_b_output_1(1:expected_age, focal_type, :, Z(:, 2) == 1, ~Z(:, 2) == 1), 3);
fe_d_2 = sum(fe_d_output_1(1:expected_age, focal_type, :, Z(:, 2) == 1, ~Z(:, 2) == 1), 3);

fe_b_2 = fe_b_2(fe_count_2 ~= 0) ./ fe_count_2(fe_count_2 ~= 0);
fe_d_2 = fe_d_2(fe_count_2 ~= 0) ./ fe_count_2(fe_count_2 ~= 0);
fe_count_2 = fe_count_2(fe_count_2 ~= 0);

fe_count_3 = sum(fe_count_output_1(1:expected_age, focal_type, :, Z(:, 3) == 1, ~Z(:, 3) == 1), 3);
fe_b_3 = sum(fe_b_output_1(1:expected_age, focal_type, :, Z(:, 3) == 1, ~Z(:, 3) == 1), 3);
fe_d_3 = sum(fe_d_output_1(1:expected_age, focal_type, :, Z(:, 3) == 1, ~Z(:, 3) == 1), 3);

fe_b_3 = fe_b_3(fe_count_3 ~= 0) ./ fe_count_3(fe_count_3 ~= 0);
fe_d_3 = fe_d_3(fe_count_3 ~= 0) ./ fe_count_3(fe_count_3 ~= 0);
fe_count_3 = fe_count_3(fe_count_3 ~= 0);

range = 5;

[counts1, edges] = histcounts(log(fe_b_1) - log(fe_d_1), -range:0.5:range);
[counts2, edges] = histcounts(log(fe_b_2) - log(fe_d_2), -range:0.5:range);
[counts3, edges] = histcounts(log(fe_b_3) - log(fe_d_3), -range:0.5:range);

center = 0.5 * (edges(1:end-1) + edges(2:end)) - 0.25;
y = [counts1; counts2; counts3]' ./ sum(sum([counts1; counts2; counts3]'));
b = bar(center, y, 1, 'stacked', 'FaceColor', 'flat', 'EdgeColor','none');

b(1).CData = hex2rgb('#03a9f4');
b(2).CData = hex2rgb('#f44336');
b(3).CData = hex2rgb('#ffc107');

b(1).FaceAlpha = 1;
b(2).FaceAlpha = 1;
b(3).FaceAlpha = 1;




ylim(y_lim)
xlim(x_lim)

if y_lab ~= ""
    ylabel(y_lab, 'FontSize', 10)
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

