function plot_error_iteration(data, colors, params)
% data structure
lambda_err = data.Y;
n_iter = data.nX;
lambda_vals = data.PLOT;
n_plot = data.nPLOT;


std_err = std(log10(lambda_err), 0, 1);
std_err = reshape(std_err, n_plot, n_iter + 1);
mean_err = mean(log10(lambda_err), 1);
mean_err = reshape(mean_err, n_plot, n_iter + 1);
upper_err = mean_err + std_err;
lower_err = mean_err - std_err;

% colors = color_profile(n_plot);
% colors = [0 0 1; 1 0 0; 0 1 0; 0 1 1; 1 0 1; 1 1 0; 0 0 0];
figure;
ax = axes;
% ax.ActivePositionProperty = 'outerposition';
for i = 1:n_plot
    fig = fill(ax, [(1:n_iter + 1)'; (n_iter + 1:-1:1)'], [10.^lower_err(i, :), 10.^rot90(upper_err(i, :), 2)],...
        colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.15);
    set(get(get(fig,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    hold on;
    string = sprintf('%.2e', lambda_vals(i)); % Convert number to a scientific notation string with 2 decimal places of precision
    stringParts = strsplit(string, 'e'); % Split the string where 'e' is
    a = str2double(stringParts(1));
    b = str2double(stringParts(2));
    if i ~= n_plot
        if a == 1 && b == 0
            lgd_txt = '$\lambda = 1$';
        elseif a == 1 && b == 1
            lgd_txt = '$\lambda = 10$';
        elseif b == 0
            lgd_txt = sprintf('$\\lambda = %.2f$', a);
        elseif a == 1
            lgd_txt = sprintf('$\\lambda = 10^{%d}$', b);
        else
            lgd_txt = sprintf('$\\lambda = %.2f \\times 10^{%d}$', a, b);
        end
    else
        if a == 1 && b == 0
            lgd_txt = 'CVX ($\lambda = 1$)';
        elseif a == 1 && b == 1
            lgd_txt = 'CVX ($\lambda = 10$)';
        elseif b == 0
            lgd_txt = sprintf('CVX ($\\lambda = %.2f$)', a);
        elseif a == 1
            lgd_txt = sprintf('CVX ($\\lambda = 10^{%d}$)', b);
        else
            lgd_txt = sprintf('CVX ($\\lambda = %.2f \\times 10^{%d}$)', a, b);
        end
    end

%     lgd_txt = sprintf('%2.0e', lambda_vals(i));
    plot(ax, 10.^mean_err(i, :), 'LineWidth', params.linewidth, 'color', colors(i, :), 'DisplayName', lgd_txt);
end

set(gca, 'YScale', 'log', 'FontSize', params.axis_fontsize, 'FontName', params.fontname);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.75 1]);
xlabel(params.xlabel_txt, 'FontName', params.fontname, 'Fontsize', params.xlabel_fontsize, 'Interpreter', params.interpreter); xlim([1 n_iter + 1]);
ylabel(params.ylabel_txt, 'Interpreter', params.interpreter, 'Fontsize', params.ylabel_fontsize, 'FontWeight', 'bold'); % ylim([0.03 3]);
lgd = legend('show', 'Location', params.legend_position);
lgd.FontSize = params.legend_fontsize;
lgd.FontName = params.fontname;
lgd.Interpreter = params.interpreter;
grid on;
grid minor;
title(params.title_txt, 'FontName', params.fontname, 'Fontsize', params.title_fontsize, 'Interpreter', params.interpreter);
axis([1 n_iter params.y_range(1) params.y_range(2)]);

% block of code to remove plot gray margin
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
