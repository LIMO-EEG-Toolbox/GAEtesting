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
    val = squeeze(mean(chan_drop(iBeta,:,:), 2));
    [grid_or_val, ~, ~]= ...
            topoplot(val,expected_chanlocs,'style','both','electrodes','off','hcolor','none','numcontour',0,'whitebk','on','noplot','on','conv','on');
    freqmap = grid_or_val(end:-1:1,:);
    freqmap(34,67)=NaN;freqmap(67,34)=NaN;freqmap(34,1)=NaN;freqmap(1,34)=NaN;
    if min(freqmap(:))<0
        freqmap = freqmap + abs(min(freqmap(:)));
    end

    figure('Color','w','NumberTitle','off','Name','limo_best_electrodes.m')
    opt = {'electrodes','on','maplimits','maxmin','verbose','off','colormap', flipud(limo_color_images(freqmap))};
    topoplot(val,expected_chanlocs,opt{:});
    title(sprintf('Channel Importance Drop Method | beta #%d', iBeta))
    exportgraphics(gcf, fullfile(plotSaveFolder, sprintf('Imp_scalp_beta_%02d_drop.png', iBeta)));
    close(gcf);

    val = squeeze(mean(chan_mean(iBeta,:,:), 2));
    [grid_or_val, ~, ~]= ...
            topoplot(val,expected_chanlocs,'style','both','electrodes','off','hcolor','none','numcontour',0,'whitebk','on','noplot','on','conv','on');
    freqmap = grid_or_val(end:-1:1,:);
    freqmap(34,67)=NaN;freqmap(67,34)=NaN;freqmap(34,1)=NaN;freqmap(1,34)=NaN;
    if min(freqmap(:))<0
        freqmap = freqmap + abs(min(freqmap(:)));
    end

    figure('Color','w','NumberTitle','off','Name','limo_best_electrodes.m')
    opt = {'electrodes','on','maplimits','maxmin','verbose','off','colormap', flipud(limo_color_images(freqmap))};
    topoplot(val,expected_chanlocs,opt{:});
    title(sprintf('Channel Importance Mean Method | beta #%d', iBeta))
    exportgraphics(gcf, fullfile(plotSaveFolder, sprintf('Imp_scalp_beta_%02d_mean.png', iBeta)));
    close(gcf);

    val = squeeze(mean(chan_zero(iBeta,:,:), 2));
    [grid_or_val, ~, ~]= ...
            topoplot(val,expected_chanlocs,'style','both','electrodes','off','hcolor','none','numcontour',0,'whitebk','on','noplot','on','conv','on');
    freqmap = grid_or_val(end:-1:1,:);
    freqmap(34,67)=NaN;freqmap(67,34)=NaN;freqmap(34,1)=NaN;freqmap(1,34)=NaN;
    if min(freqmap(:))<0
        freqmap = freqmap + abs(min(freqmap(:)));
    end

    figure('Color','w','NumberTitle','off','Name','limo_best_electrodes.m')
    opt = {'electrodes','on','maplimits','maxmin','verbose','off','colormap', flipud(limo_color_images(freqmap))};
    topoplot(val,expected_chanlocs,opt{:});
    title(sprintf('Channel Importance Zero Method | beta #%d', iBeta))
    exportgraphics(gcf, fullfile(plotSaveFolder, sprintf('Imp_scalp_beta_%02d_zero.png', iBeta)));
    close(gcf);
end