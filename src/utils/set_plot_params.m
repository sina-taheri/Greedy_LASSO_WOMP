function params = set_plot_params(i_experiment, title_txt, legend_txt)

params.fontname = 'Times';
params.interpreter = 'latex';

params.title_txt = title_txt;
params.title_fontsize = 36;

params.marker_size = 8;
params.linewidth = 2;

params.axis_fontsize = 30;

params.xlabel_fontsize = 36;
params.ylabel_fontsize = 36;

switch i_experiment
    case 1
        params.xlabel_txt = '$\lambda$';
        params.ylabel_txt = 'Relative $\ell^2$-error';

        params.legend_fontsize = 36;
        params.legend_position = 'northwest';
        params.legend_txt = legend_txt;

        params.y_range = [10^-4.5 10^1];
    case 2
        params.xlabel_txt = '$\lambda$';
        params.ylabel_txt = 'Relative $L_\varrho^2(\mathcal{D})$-error';

        params.legend_fontsize = 36;
        params.legend_position = 'northwest';
        params.legend_txt = legend_txt;

        params.y_range = [10^-4 10^1];
    case 3
        params.xlabel_txt = '$\lambda$';
        params.ylabel_txt = 'Relative $\ell^2$-error';
        
        params.legend_fontsize = 36;
        params.legend_position = 'southeast';
        params.legend_txt = legend_txt;

        params.y_range = [10^-4.5 10^1];
    case 4
        params.xlabel_txt = 'Iteration Number';
        params.ylabel_txt = 'Relative $\ell^2$-error';

        params.legend_fontsize = 32;
        params.legend_position = 'northwest';

        params.y_range = [10^-4.5 10^2];
    case 5
        params.xlabel_txt = 'Iteration Number';
        params.ylabel_txt = 'Relative $ L_\varrho^2(\mathcal{D})\ $-error ';

        params.legend_fontsize = 32;
        params.legend_position = 'northwest';

        params.y_range = [10^-3 10^0.5];
end
        