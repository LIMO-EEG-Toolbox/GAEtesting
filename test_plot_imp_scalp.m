% Set path
curr_path = pwd();
analysis_type = 'N170';
chanlocs_path = fullfile(curr_path, 'limo_chanlocs.mat');
maskingRes_path = fullfile(curr_path, analysis_type, 'Masking', 'Channel Masking');
chan_drop_path = fullfile(maskingRes_path, 'imp_chan_drop.mat');
chan_mean_path = fullfile(maskingRes_path, 'imp_chan_mean.mat');
chan_zero_path = fullfile(maskingRes_path, 'imp_chan_zero.mat');

% Create plot save folder
plotSaveFolder = fullfile(maskingRes_path, 'Scalp Plots');
if ~exist(plotSaveFolder, 'dir')
    mkdir(plotSaveFolder);
end

% Load the adjacency matrix and data
chan_drop = load(chan_drop_path).all_imp_chan_drop;
chan_mean = load(chan_mean_path).all_imp_chan_mean;
chan_zero = load(chan_zero_path).all_imp_chan_zero;
expected_chanlocs = load(chanlocs_path).expected_chanlocs;

% Generate scalp plots for each beta
[nBeta, nSubject, nChan] = size(chan_drop);

for iBeta = 1:nBeta
    val_drop = squeeze(mean(chan_drop(iBeta,:,:), 2));
    title_drop = sprintf('Channel Importance Drop Method | beta #%d', iBeta);
    filename_drop = sprintf('Imp_scalp_beta_%02d_drop.png', iBeta);
    plot_chan_importance(val_drop, expected_chanlocs, ...
                         title_drop, plotSaveFolder, filename_drop);

    val_mean = squeeze(mean(chan_mean(iBeta,:,:), 2));
    title_mean = sprintf('Channel Importance Mean Method | beta #%d', iBeta);
    filename_mean = sprintf('Imp_scalp_beta_%02d_mean.png', iBeta);
    plot_chan_importance(val_mean, expected_chanlocs, ...
                         title_mean, plotSaveFolder, filename_mean);

    val_zero = squeeze(mean(chan_zero(iBeta,:,:), 2));
    title_zero = sprintf('Channel Importance Zero Method | beta #%d', iBeta);
    filename_zero = sprintf('Imp_scalp_beta_%02d_zero.png', iBeta);
    plot_chan_importance(val_zero, expected_chanlocs, ...
                         title_zero, plotSaveFolder, filename_zero);
end

function plot_chan_importance(val, expected_chanlocs, ...
                         title_str, plotSaveFolder, filename)

    val_centered = val - mean(val);
    c = max(abs(val_centered));
    figure('Color','w','NumberTitle','off','Name','limo_best_electrodes.m');
    opt = {'electrodes','on', ...
           'maplimits',[-c c], ...
           'verbose','off', ...
           'colormap', flipud(limo_color_images(val_centered))};

    topoplot(val_centered, expected_chanlocs, opt{:});
    %colorbar;
    title(title_str);
    exportgraphics(gcf, fullfile(plotSaveFolder, filename));
    close(gcf);
end