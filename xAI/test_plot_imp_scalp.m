% Set path
curr_path = pwd();
analysis_type = 'P3';
chanlocs_path = fullfile(curr_path, 'datain', 'limo_chanlocs.mat');
maskingRes_path = fullfile(curr_path, analysis_type, 'Ablation', 'Channel Ablation');
chan_mean_path = fullfile(maskingRes_path, 'imp_chan_drop.mat');

% Create plot save folder
plotSaveFolder = fullfile(maskingRes_path, 'Scalp Plots');
if ~exist(plotSaveFolder, 'dir')
    mkdir(plotSaveFolder);
end

% Load the adjacency matrix and data
chan_mean = load(chan_mean_path).all_imp_chan_drop;
expected_chanlocs = load(chanlocs_path).expected_chanlocs;

% Generate scalp plots for each beta
[nBeta, nSubject, nChan] = size(chan_mean);

chan_mean_norm = zeros(size(chan_mean));
eps_val = 1e-12;

for iBeta = 1:nBeta
    X = squeeze(chan_mean(iBeta,:,:));
    denom = max(X(:));
    denom = max(denom, eps_val);
    chan_mean_norm(iBeta,:,:) = X ./ denom;
end

beta_chan_mean = squeeze(mean(chan_mean_norm, 2));
final_chan_imp = mean(beta_chan_mean, 1)';
filename  = sprintf('%s_Imp_scalp.png', analysis_type);
title_str = sprintf('%s | Aggregated Channel Importance', analysis_type);
fprintf('Saved scalp plot: %s\n', fullfile(plotSaveFolder, filename));

plot_chan_importance(final_chan_imp, expected_chanlocs, title_str, plotSaveFolder, filename);

function plot_chan_importance(val, expected_chanlocs, ...
                         title_str, plotSaveFolder, filename)

    val_centered = val - mean(val);
    c = max(abs(val_centered));
    figure('Color','w','NumberTitle','off','Name','limo_best_electrodes.m');
    val_plot = val;
    cmax = max(val_plot);
    opt = {'electrodes','on', ...
       'maplimits',[0 cmax], ...
       'verbose','off', ...
       'colormap', limo_color_images(val_plot)};

    topoplot(val_centered, expected_chanlocs, opt{:});
    %colorbar;
    title(title_str);
    exportgraphics(gcf, fullfile(plotSaveFolder, filename));
    close(gcf);
end