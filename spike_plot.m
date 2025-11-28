% map_cluster_properties_spatial.m
% Visualize per-cluster properties on the probe at each cluster's BEST CHANNEL.
% Also shows all metrics for a selected cluster and its template waveform.
% clc;close all;clear
% HOW TO USE:
  cfg = struct;
  cfg.phyDir = 'Z:/HarveyLab/Tier1/Kevin/Videos/KM49-50/251120_g0/Spike_Sorting';
  cfg.propertyName = 'amp';        % any column from cluster_info.tsv / cluster_*.tsv / computed n_spikes
  cfg.agg = 'max';                 % how to combine if multiple clusters share a channel: 'max'|'mean'|'count'|'first'
  cfg.onlyIds = [];                % restrict to specific cluster IDs ([])
  cfg.excludeIds = [];             % exclude some cluster IDs ([])
  cfg.saveFig = false;             % save overview figure
  cfg.outPng = 'property_map.png';
  % [templates, Tmetrics]=map_cluster_properties_spatial(cfg);
%
% Inside the UI: click a point (cluster) to update the right-hand panels (metrics & waveform).
%
% REQUIREMENTS:
%   - readNPY.m on MATLAB path
%   - Kilosort/Phy files: templates.npy, spike_clusters.npy, spike_templates.npy
%     optional: templates_ind.npy, channel_positions.npy (or chanMap.mat), channel_ids.npy/channel_map.npy
%     optional metrics: cluster_info.tsv and any cluster_*.tsv (with 'cluster_id' column)

% function [templates, Tmetrics]=map_cluster_properties_spatial(cfg)
assert(exist('readNPY','file')==2, 'readNPY.m not found on MATLAB path.');

% ------- defaults -------
need = @(f) isfield(cfg,f) && ~isempty(cfg.(f));
if ~need('phyDir'), error('Set cfg.phyDir'); end
if ~need('propertyName'), cfg.propertyName = ''; end
if ~need('agg'), cfg.agg = 'max'; end
if ~need('onlyIds'), cfg.onlyIds = []; end
if ~need('excludeIds'), cfg.excludeIds = []; end
if ~need('saveFig'), cfg.saveFig = false; end
if ~need('outPng'), cfg.outPng = 'property_map.png'; end

phy = fullfile(cfg.phyDir,'phy');

% ------- load arrays -------
templates_raw      = readNPY(fullfile(phy,'templates.npy'));          % (T,time,chan) or (T,chan,time) or (T,time)
spike_clusters = readNPY(fullfile(phy,'spike_clusters.npy'));     % (Nspikes,)
spike_templates= readNPY(fullfile(phy,'spike_templates.npy'));    % (Nspikes,)
templates_ind  = tryRead(fullfile(phy,'template_ind.npy'));      % (T,chanSubset) global channel ids
channel_ids    = firstHitNPY(phy, {'channel_ids.npy','channel_map.npy'});  % optional
spike_times    = readNPY(fullfile(phy,'spike_times.npy'));
amplitudes    = readNPY(fullfile(phy,'amplitudes.npy'));
chanPos        = loadChannelPositions(phy);                       % (Nchan,2) x,y (µm)

% orient templates → (nTemplates, nTime, nChan)
templates = orientTemplates(templates_raw);

% metrics table (cluster_info.tsv + cluster_*.tsv) + computed n_spikes
% Tmetrics = buildMetricsTable(phy, spike_clusters);
T = parquetread(fullfile(cfg.phyDir,'quality_metrics.parquet'));   % table of metrics
% crt="Tmetrics.snr > 2.0 & Tmetrics.fr > 0.05 & Tmetrics.nn_hit_rate > 0.5 & Tmetrics.isi_violations_ratio<1 & Tmetrics.amplitude_cutoff<0.1 & Tmetrics.presence_ratio>0.9";
crt="T.snr>2.0 & T.firing_rate>0.01 & T.isi_violations_ratio<1 & T.presence_ratio>0.9 & T.sync_spike_2<0.15"
allIds = find(eval(crt));

% filter cluster ids (based on presence in spike_clusters, plus user lists)
% allIds = unique(double(spike_clusters(:)));
if ~isempty(cfg.onlyIds), allIds = intersect(allIds, cfg.onlyIds(:)); end
if ~isempty(cfg.excludeIds), allIds = setdiff(allIds, cfg.excludeIds(:)); end

% compute best channel + representative template per cluster
C = computeClusterSummaries(allIds, templates, templates_ind, spike_clusters, spike_templates, channel_ids, chanPos);

% attach metrics row to C
% C = attachMetrics(C, Tmetrics);
%
norm_template=normalize(templates(allIds,:,1),2,"range");
xy={C.best_xy};
arr_xy=cell2mat(xy');
figure;
cm=colormap("lines");
for i=1: length(C)
    hold on
    
    
plot([1:120]/4 + arr_xy(i,1),norm_template(i,:)*20 + arr_xy(i,2),"Color",cm(ceil(allIds(i)/length(templates_ind)*256),:),"LineWidth",1)
end
scatter(chanPos(:,1)+20,chanPos(:,2),[],'k')
title(strcat(crt,"  n=",string(length(allIds))),"Interpreter","none")
% ylim([0,2100])
xlim([0 900])



% pick property column
propName = cfg.propertyName;
if isempty(propName)
    % choose a reasonable default
    if any(strcmpi('snr', Tmetrics.Properties.VariableNames))
        propName = 'snr';
    elseif any(strcmpi('amplitude', Tmetrics.Properties.VariableNames))
        propName = 'amplitude';
    else
        propName = 'n_spikes';
    end
end
assert(any(strcmp(propName, Tmetrics.Properties.VariableNames)), ...
    'Property "%s" not found. Available columns:\n%s', propName, strjoin(Tmetrics.Properties.VariableNames, ', '));

% aggregate clusters to channels for the overview map
[chXY, chVal, chClusters] = aggregateToChannels(C, propName, cfg.agg);

% ------- figure layout -------
fig = figure('Name','Cluster property map (best channels)','Color','w','Position',[80 80 1200 600]);
t = tiledlayout(fig, 2, 3, 'TileSpacing','compact','Padding','compact');
title(t, sprintf('Property: %s   |   Agg: %s', propName, cfg.agg), 'FontWeight','bold');

% (1) probe map with property color (all clusters aggregated to channels)
ax1 = nexttile(t, [2 2]);
hold(ax1,'on');
% draw all physical channels in light gray
if ~isempty(chanPos)
    scatter(ax1, chanPos(:,1), chanPos(:,2), 10, [0.85 0.85 0.85], 'filled', 'MarkerEdgeColor','none');
end
% draw channel-level aggregated property
if ~isempty(chXY)
    s = scatter(ax1, chXY(:,1), chXY(:,2), 60, chVal, 'filled', 'MarkerEdgeColor','k');
    colormap(ax1, parula); colorbar(ax1);
    caxis(ax1, robustLimits(chVal));
end
axis(ax1,'equal'); grid(ax1,'on'); box(ax1,'on');
set(ax1,'YDir','reverse'); % common convention for probe plots
xlabel(ax1,'x (µm)'); ylabel(ax1,'y (µm)');
title(ax1, 'Best-channel property (channel-aggregated)');

% interactive picking: show one cluster's details when clicking near a channel
set(fig, 'WindowButtonDownFcn', @(~,~) onClick(ax1, C, chXY, chClusters));

% (2) right top: waveform for selected cluster
ax2 = nexttile(t, 1);
title(ax2, 'Template waveform @ best channel'); grid(ax2,'on'); box(ax2,'on');

% (3) right bottom: metrics for selected cluster
ax3 = nexttile(t, 1);
axis(ax3,'off'); title(ax3,'Cluster metrics');

% Store handles & data for callbacks
S.ax1 = ax1; S.ax2 = ax2; S.ax3 = ax3; S.C = C; S.chanPos = chanPos;
guidata(fig, S);

% initial selection: pick cluster with max property
[~,imax] = max(chVal);
if ~isempty(imax) && ~isempty(chClusters)
    cid0 = chClusters{imax}(1);
    showClusterDetails(fig, cid0);
end

if cfg.saveFig
    exportgraphics(fig, cfg.outPng, 'Resolution', 200);
end

% end % main

%% ============================ helpers ============================

function arr = tryRead(p)
    if exist(p,'file'), arr = readNPY(p); else, arr = []; end
end

function arr = firstHitNPY(root, names)
    arr = [];
    for k = 1:numel(names)
        p = fullfile(root, names{k});
        if exist(p,'file'), arr = readNPY(p); return; end
    end
end

function chanPos = loadChannelPositions(root)
    chanPos = [];
    p = fullfile(root,'channel_positions.npy');
    if exist(p,'file')
        chanPos = readNPY(p);
        return
    end
    % try Kilosort chanMap.mat
    pm = fullfile(root,'chanMap.mat');
    if exist(pm,'file')
        S = load(pm);
        % common fields: xcoords, ycoords
        if isfield(S,'xcoords') && isfield(S,'ycoords')
            chanPos = [S.xcoords(:), S.ycoords(:)];
            return
        end
    end
end

function templates = orientTemplates(templates)
% Return (nTemplates, nTime, nChan) from (T,time,chan) | (T,chan,time) | (T,time)
    sz = size(templates);
    if ndims(templates)==3
        T = sz(1); d2 = sz(2); d3 = sz(3);
        if d2 <= 512
            % assume (T, time, chan)
        else
            % assume (T, chan, time)
            templates = permute(templates,[1 3 2]);
        end
    elseif ismatrix(templates)
        templates = reshape(templates, sz(1), sz(2), 1);
    else
        error('Unexpected templates ndim: %d', ndims(templates));
    end
end

function [bestGlobalCh, bestLocal] = bestChannelByPTP(T_time_by_chan, chInds)
    if isempty(chInds), chInds = 0:(size(T_time_by_chan,2)-1); end
    ptp = max(T_time_by_chan,[],1) - min(T_time_by_chan,[],1);
    [~,bestLocal] = max(ptp);
    bestGlobalCh = chInds(bestLocal);
end

function repT = representativeTemplateForCluster(spike_templates, spikeIdx)
    st = double(spike_templates(spikeIdx));
    u = unique(st);
    [~,~,ic] = unique(st);
    cnt = accumarray(ic,1);
    [~,imax] = max(cnt);
    repT = u(imax); % 0-based
end

function Tmetrics = buildMetricsTable(phy, spike_clusters)
    sc = double(spike_clusters(:));
    [u,~,ic] = unique(sc);
    cnt = accumarray(ic,1);
    Tcount = table(u, cnt, 'VariableNames', {'cluster_id','n_spikes'});

    Tmetrics = Tcount;
    pinfo = fullfile(phy,'cluster_info.tsv');
    if exist(pinfo,'file')
        try
            Tin = readtable(pinfo,'FileType','text','Delimiter','\t');
            if ~ismember('cluster_id', Tin.Properties.VariableNames)
                if ismember('id', Tin.Properties.VariableNames)
                    Tin.Properties.VariableNames{strcmp(Tin.Properties.VariableNames,'id')} = 'cluster_id';
                end
            end
            Tmetrics = outerjoin(Tmetrics, Tin, 'Keys','cluster_id', 'MergeKeys',true);
        catch ME
            warning('cluster_info.tsv read failed: %s', ME.message);
        end
    end
    % d = dir(fullfile(phy,'cluster_*.tsv'));
    % for k = 1:numel(d)
    %     if strcmp(d(k).name,'cluster_info.tsv'), continue; end
    %     try
    %         T = readtable(fullfile(phy,d(k).name),'FileType','text','Delimiter','\t');
    %         if ismember('cluster_id', T.Properties.VariableNames)
    %             Tmetrics = outerjoin(Tmetrics, T, 'Keys','cluster_id', 'MergeKeys',true);
    %         end
    %     catch ME
    %         warning('%s read failed: %s', d(k).name, ME.message);
    %     end
    % end
    Tmetrics = sortrows(Tmetrics,'cluster_id');
end

function C = computeClusterSummaries(clusterIds, templates, templates_ind, spike_clusters, spike_templates, channel_ids, chanPos)
    C = struct('cluster_id',{},'rep_template',{},'best_global_ch',{},'best_local_ch',{},'best_xy',{},'ptp',{},'row',[]);
    for i = 1:numel(clusterIds)
        cid = clusterIds(i);
        spikeIdx = find(double(spike_clusters)==cid);
        if isempty(spikeIdx), continue; end
        repT = representativeTemplateForCluster(spike_templates, spikeIdx);
        T = squeeze(templates(repT+1,:,:)); % (time, chansub)
        if isvector(T), T = T(:); end
        if size(T,1) < size(T,2), T = T.'; end

        if ~isempty(templates_ind)
            chInds = templates_ind(repT+1,:);
            chInds = chInds(:).';
            nc = size(T,2);
            if numel(chInds) < nc
                chInds = [chInds, 0:(nc-numel(chInds)-1)];
            elseif numel(chInds) > nc
                chInds = chInds(1:nc);
            end
        else
            chInds = 0:(size(T,2)-1);
        end

        [bestCh, bestLoc] = bestChannelByPTP(T, chInds);
        ptp = max(T(:,bestLoc)) - min(T(:,bestLoc));
        xy = [NaN NaN];
        if ~isempty(chanPos) && (bestCh+1) <= size(chanPos,1), xy = chanPos(bestCh+1,:); end
        if ~isempty(channel_ids) && (bestCh+1) <= numel(channel_ids)
            % optional: remap bestCh label to provided channel_ids value
        end

        C(end+1).cluster_id = cid; %#ok<AGROW>
        C(end).rep_template = repT;
        C(end).best_global_ch = bestCh;
        C(end).best_local_ch  = bestLoc;
        C(end).best_xy = xy;
        C(end).ptp = ptp;
    end
end

function C = attachMetrics(C, Tmetrics)
    if isempty(C), return; end
    % attach the entire row (all properties) from Tmetrics to C(i).row
    m = containers.Map('KeyType','double','ValueType','double');
    ids = Tmetrics.cluster_id;
    for i = 1:numel(C)
        cid = C(i).cluster_id;
        idx = find(ids==cid, 1, 'first');
        if ~isempty(idx)
            C(i).row = Tmetrics(idx,:);
        else
            C(i).row = Tmetrics(1,:); C(i).row(1:end,:) = {[]}; % empty row
        end
    end
end

function [XY, V, chClusters] = aggregateToChannels(C, propName, agg)
    % Build per-channel value by aggregating chosen property across clusters that share the channel
    XY = []; V = []; chClusters = {};
    if isempty(C), return; end

    % pull values
    chanMap = containers.Map('KeyType','double','ValueType','any');
    for i = 1:numel(C)
        cid = C(i).cluster_id;
        if isempty(C(i).row)
            val = NaN;
        else
            row = C(i).row;
            if ismember(propName, row.Properties.VariableNames)
                val = row.(propName);
                if istable(val), val = val{1}; end
            else
                val = NaN;
            end
        end
        ch = double(C(i).best_global_ch);
        if ~isKey(chanMap, ch)
            chanMap(ch) = struct('vals', val, 'xy', C(i).best_xy, 'cids', cid);
        else
            s = chanMap(ch);
            s.vals = [s.vals, val]; %#ok<AGROW>
            s.cids = [s.cids, cid]; %#ok<AGROW>
            chanMap(ch) = s;
        end
    end

    % aggregate
    ks = sort(cell2mat(keys(chanMap)));
    XY = zeros(numel(ks),2); V = nan(numel(ks),1); chClusters = cell(numel(ks),1);
    for i = 1:numel(ks)
        s = chanMap(ks(i));
        XY(i,:) = s.xy;
        chClusters{i} = s.cids;
        vals = s.vals;
        switch lower(agg)
            case 'mean', V(i) = mean(vals(~isnan(vals)));
            case 'max',  V(i) = max(vals(~isnan(vals)));
            case 'count',V(i) = numel(vals);
            case 'first',V(i) = vals(find(~isnan(vals),1,'first'));
            otherwise,   V(i) = max(vals(~isnan(vals)));
        end
        if isempty(V(i)) || isnan(V(i)), V(i) = NaN; end
    end
end

function lim = robustLimits(x)
    x = x(:); x = x(~isnan(x) & isfinite(x));
    if isempty(x), lim = [0 1]; return; end
    p = prctile(x,[2 98]);
    if p(1)==p(2), p = [min(x) max(x)]; end
    lim = p;
end

function onClick(ax, C, chXY, chClusters)
    S = guidata(ax.Parent);
    if isempty(chXY), return; end
    cp = get(ax, 'CurrentPoint'); p = cp(1,1:2);
    % nearest channel dot
    d2 = sum((chXY - p).^2,2);
    [~,i] = min(d2);
    if isempty(chClusters) || i<1 || i>numel(chClusters) || isempty(chClusters{i}), return; end
    % if multiple clusters share channel, pick the first; you can change this logic
    cid = chClusters{i}(1);
    showClusterDetails(ax.Parent, cid);
end

function showClusterDetails(fig, cid)
    S = guidata(fig);
    C = S.C;
    i = find([C.cluster_id]==cid,1,'first'); if isempty(i), return; end

    % --- waveform panel ---
    ax2 = S.ax2; cla(ax2);
    repT = C(i).rep_template;
    % NOTE: we don't have templates here; re-read minimal slice from guidata (store if needed)
    % For simplicity, draw a schematic from C(i).ptp only if templates not cached.
    % Better: cache templates in S; to keep memory sane with big arrays we plot from disk once.
    % Here we regenerate from disk quickly:
    phy = fileparts(fileparts(mfilename('fullpath'))); %#ok<NASGU>
    % Instead of reloading, store a function handle in S to get waveform:
    if ~isfield(S,'getWaveformFcn')
        warning('Waveform accessor not found; attach it if you want live waveforms here.');
        text(ax2, .5,.5, sprintf('Cluster %d\nPTP=%.3g (no live template cached)', cid, C(i).ptp), ...
            'HorizontalAlignment','center'); axis(ax2,'off');
    else
        [tvec, wf, label] = S.getWaveformFcn(C(i));
        plot(ax2, tvec, wf, 'LineWidth', 1.25);
        title(ax2, sprintf('Cluster %d — Template @ best ch %s', cid, label));
        xlabel(ax2,'Samples'); ylabel(ax2,'Amplitude (a.u.)'); grid(ax2,'on'); box(ax2,'on');
    end

    % --- metrics panel ---
    ax3 = S.ax3; cla(ax3); axis(ax3,'off');
    row = C(i).row;
    y = 0.95; dy = 0.045;
    text(ax3, 0.02, y, sprintf('cluster_id: %d', cid), 'FontWeight','bold','Interpreter','none'); y = y - dy;
    if ~isempty(row)
        vn = row.Properties.VariableNames;
        for k = 1:numel(vn)
            v = vn{k};
            if strcmp(v,'cluster_id'), continue; end
            val = row.(v);
            if istable(val), val = val{1}; end
            if iscell(val), val = val{1}; end
            if isnumeric(val)
                s = sprintf('%s: %g', v, val);
            elseif isstring(val) || ischar(val)
                s = sprintf('%s: %s', v, string(val));
            else
                s = sprintf('%s: [..]', v);
            end
            text(ax3, 0.02, y, s, 'Interpreter','none');
            y = y - dy; if y < 0.05, break; end
        end
    end
    drawnow;
end

%% After creating the figure, attach a waveform accessor:
fig = gcf; S = guidata(fig);

% Cache what we need to read waveforms quickly
S.templates = readNPY(fullfile(cfg.phyDir,'templates.npy'));
S.templates  = orientTemplates(S.templates);    % reuse the helper above
S.templates_ind = tryRead(fullfile(cfg.phyDir,'templates_ind.npy'));
S.channel_ids   = firstHitNPY(cfg.phyDir, {'channel_ids.npy','channel_map.npy'});

S.getWaveformFcn = @(Ci) localGetWF(Ci, S.templates, S.templates_ind, S.channel_ids);
guidata(fig, S);

function [tvec, wf, label] = localGetWF(Ci, templates, templates_ind, channel_ids)
    T = squeeze(templates(Ci.rep_template+1,:,:));
    if isvector(T), T = T(:); end
    if size(T,1) < size(T,2), T = T.'; end
    if ~isempty(templates_ind)
        chInds = templates_ind(Ci.rep_template+1,:); chInds = chInds(:).';
        if numel(chInds) < size(T,2), chInds = [chInds, 0:(size(T,2)-numel(chInds)-1)]; end
        if numel(chInds) > size(T,2), chInds = chInds(1:size(T,2)); end
    else
        chInds = 0:(size(T,2)-1);
    end
    ptp = max(T,[],1)-min(T,[],1);
    [~,k] = max(ptp);
    wf = T(:,k);
    tvec = (0:numel(wf)-1).';
    bestCh = chInds(k);
    if ~isempty(channel_ids) && (bestCh+1)<=numel(channel_ids)
        label = num2str(channel_ids(bestCh+1));
    else
        label = num2str(bestCh);
    end
end
