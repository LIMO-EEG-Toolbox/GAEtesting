function run_limo_ttest(work_dir, limo_tools_dir, chanlocs_mat)
% Cluster Python two-sample t-maps with LIMO (t^2, min 2 channels).
% Does not call limo_random_robust (needs Statistics Toolbox tcdf).

if nargin < 3
    error('usage: run_limo_ttest(work_dir, limo_tools_dir, chanlocs_mat)')
end

addpath(limo_tools_dir);
addpath(fullfile(limo_tools_dir, 'limo_cluster_functions'));
addpath(fullfile(limo_tools_dir, 'external'));
addpath(fullfile(limo_tools_dir, 'external', 'psom'));
% IPT / SPM mex missing on this Mac; limo_sims/bwlabeln.m is 2D 4-connected.

work_dir = char(work_dir);
S = load(fullfile(work_dir, 'bundle.mat'));
loc = load(char(chanlocs_mat));

LIMO = struct();
LIMO.Level = 2;
LIMO.Analysis = 'Time';
LIMO.Type = 'Channels';
LIMO.data.chanlocs = loc.expected_chanlocs;
LIMO.data.neighbouring_matrix = loc.channeighbstructmat;
LIMO.dir = work_dir;

names = {'mean', 'yuen', 'weighted'};
summary = struct();
ch = double(S.planted_ch(:));
ch1 = [];
ch2 = [];
if isfield(S, 'planted_ch1')
    ch1 = double(S.planted_ch1(:));
end
if isfield(S, 'planted_ch2')
    ch2 = double(S.planted_ch2(:));
end
for i = 1:numel(names)
    n = names{i};
    t = S.([n '_t']);
    p = S.([n '_p']);
    bootT = S.([n '_bootT']);
    bootP = S.([n '_bootP']);
    [mask, cluster_pval, max_th] = limo_clustering( ...
        t.^2, p, bootT.^2, bootP, LIMO, 2, 0.05, 0);
    summary.(n) = pack_cluster(mask, cluster_pval, max_th, ch, ch1, ch2);
end

save(fullfile(work_dir, 'limo_cluster_summary.mat'), 'summary');
fid = fopen(fullfile(work_dir, 'limo_cluster_summary.json'), 'w');
fwrite(fid, jsonencode(summary));
fclose(fid);
end


function s = pack_cluster(mask, cluster_pval, max_th, ch, ch1, ch2)
if nargin < 5
    ch1 = [];
end
if nargin < 6
    ch2 = [];
end
if isempty(mask)
    mask = zeros(0);
end
s.n_sig_cells = nnz(mask);
s.max_th = max_th;
s.planted_any_cluster = false;
s.planted_n_sig = 0;
s.n1_any = false;
s.n1_n_sig = 0;
s.n2_any = false;
s.n2_n_sig = 0;
if ~isempty(mask) && ndims(mask) >= 2
    if ~isempty(ch)
        s.planted_any_cluster = any(mask(ch, :), 'all');
        s.planted_n_sig = nnz(mask(ch, :));
    end
    if ~isempty(ch1)
        s.n1_any = any(mask(ch1, :), 'all');
        s.n1_n_sig = nnz(mask(ch1, :));
    end
    if ~isempty(ch2)
        s.n2_any = any(mask(ch2, :), 'all');
        s.n2_n_sig = nnz(mask(ch2, :));
    end
end
if isempty(cluster_pval)
    s.min_cluster_p = NaN;
else
    s.min_cluster_p = min(cluster_pval(:));
end
end
