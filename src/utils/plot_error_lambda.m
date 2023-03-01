function plot_error_lambda(data, color_profile, params)
% data structure
Y = data.Y;
X = data.X;
nX = data.nX;
PLOT = data.PLOT;
nPLOT = data.nPLOT;

figure;
ax = axes;
for i = 1:nPLOT
    hbp = boxplot(Y{i}, X, 'color', color_profile(i, :));
    h = findobj(hbp, 'tag', 'Outliers');
    for iH = 1:length(h)
        h(iH).MarkerEdgeColor = color_profile(i, :);
    end
    hAx = gca;
    xtk = hAx.XTick; 
    hold on;
    plot(xtk, median(Y{i}, 1), 'Marker', '*', 'MarkerSize', params.marker_size, 'color', color_profile(i, :), 'LineWidth', params.linewidth, 'DisplayName', params.legend_txt{i});
    hold on;
end

str = sprintf('%.2e', X);
str_mat = reshape(str, [8 nX])';
num_parts = str2double(split(string(str_mat), 'e'));
xtk = cell(1, nX);
for i = 1:nX
    a = num_parts(i, 1);
    if a == 1
        b = num_parts(i, 2);
        xtk{i} = sprintf('$10^{%d}$', b);
    end
end
ax.XTickLabels = xtk;
% ax.YTickMode = 'manual';
ax.TickLabelInterpreter = params.interpreter;
set(gca, 'YScale', 'log', 'FontSize', params.axis_fontsize, 'FontName', params.fontname);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.75 1]);
xlabel(params.xlabel_txt, 'Fontsize', params.xlabel_fontsize, 'Interpreter', params.interpreter);
ylabel(params.ylabel_txt, 'Fontsize', params.ylabel_fontsize, 'Interpreter', params.interpreter);
lgd = legend('show', 'Location', params.legend_position);
lgd.FontSize = params.legend_fontsize;
lgd.FontName = params.fontname;
lgd.Interpreter = params.interpreter;
grid on;
grid minor;
title(params.title_txt, 'FontName', params.fontname, 'Fontsize', params.title_fontsize, 'Interpreter', params.interpreter);
axis([1 nX params.y_range(1) params.y_range(2)]);

% block of code to remove plot gray margin
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
