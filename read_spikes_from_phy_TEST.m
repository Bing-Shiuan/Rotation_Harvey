%% This is a test script to read sorted units from phy

% The info I want is:
% - spike times
% - spike group
% - some quality metrics?
% - where the neuron is on the probe so I can localize in the brain
% all of this is in cluster_info.tsv

% however, on this test, I screwed up and did not save the quality metrics
% only, so I'll have to quality filter here


% set paths
cd('D:\ANALYSIS\Kevin')
addpath('D:\ANALYSIS\Kevin\Utilities')

% set folders
phy_folder = 'D:\DATA\Kevin\KM33\241219_g0\Spike_Sorting\phy';
TTL_folder = 'D:\DATA\Kevin\KM33\241219_g0\TTLs\TTLs.mat';
beh_folder = 'D:\DATA\Kevin\KM33\241219*.dat';

% from cluster
phy_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM33\250113_g0\Spike_Sorting\phy';
TTL_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM33\250113_g0\TTLs\TTLs.mat';
beh_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM33\250113*.dat';


phy_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM37\250407_g0\Spike_Sorting\phy';
TTL_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM37\250407_g0\TTLs\TTLs.mat';
beh_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM37\250407*.dat';
mouse_name = 'KM37';

% new km 41?
phy_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM41\250416_g0\Spike_Sorting\phy';
TTL_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM41\250416_g0\TTLs\TTLs.mat';
beh_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM41\250416*.dat';
mouse_name = 'KM41';

% KM 42 VCA1
phy_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM42\250727_g0\Spike_Sorting\phy';
TTL_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM42\250727_g0\TTLs\TTLs.mat';
beh_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM42\250727*.dat';
mouse_name = 'KM42';

% KM 43 in BLA + HC, though in reality I totally missed and am in a fiber
% tract
phy_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM43\250727_g0\Spike_Sorting\phy';
TTL_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM43\250727_g0\TTLs\TTLs.mat';
beh_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM43\250727*.dat';
mouse_name = 'KM43';

% KM47 since implant crap, wonder if will get better?
phy_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM47\250815_g0\Spike_Sorting\phy';
TTL_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM47\250815_g0\TTLs\TTLs.mat';
beh_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM47\250815*.dat';
mouse_name = 'KM47';

%%
% set paths
clc; clear; close all;
addpath(genpath('Z:\HarveyLab\Tier1\Kevin\Analysis\20250718_backup_Cindys_PC\Utilities'))


phy_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM41\250416_g0\Spike_Sorting\phy';
TTL_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM41\250416_g0\TTLs\TTLs.mat';
beh_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM41\250416*.dat';
mouse_name = 'KM33';
%% Load data

% load unit info
spike_times = readNPY(fullfile(phy_folder, 'spike_times.npy'));
spike_clusters = readNPY(fullfile(phy_folder, 'spike_clusters.npy'));
% spike times are in samples!!!
% - so to convert to time (s), divide by sampling frequency.
 % 92,141,563 samples 


% load auxiliary info (depth, quality metrics)
try % if error occur, run catch
    T = readtable(fullfile(phy_folder, 'cluster_info.tsv'),...
        'FileType','text','Delimiter','\t');
    % this only works if I run phy
    % - better to save from my original spike sorting code...
catch
    disp('loading quality metrics from parquet?')
    qc_path = fullfile(phy_folder(1:end-4), 'quality_metrics.parquet');
    T = parquetread(qc_path);
    % um
    T = renamevars(T, 'firing_rate', 'fr');
    % need to get depth also
    channel_positions = readNPY(fullfile(phy_folder, 'channel_positions.npy'));
    channel_groups = readNPY(fullfile(phy_folder, 'channel_groups.npy'));
    channel_map = readNPY(fullfile(phy_folder, 'channel_map.npy'));
    primary_channels = readNPY(fullfile(phy_folder(1:end-4), 'primary_channels.npy'));
    % map primary channels
    % primary channels is a 1xunit array of which NP channel
    % channel map is a 1xchannel array
    % channel_positions is from the probe. the first colum is xcoordinates,
    % the second column is y coordinate (in microns)
    depth = channel_positions(primary_channels+1 , 2); % want y coordinte
    
end

%T.depth;
%T.cluster_id

% load behavioral data
beh_path = dir(beh_folder);
info = createBEHstruct_nonsocial(mouse_name, beh_path);

%% double check if depth
% does T.depth (from phy) match what I figured out from saved data?

% need to get depth also
channel_positions = readNPY(fullfile(phy_folder, 'channel_positions.npy'));
primary_channels = readNPY(fullfile(phy_folder(1:end-4), 'primary_channels.npy'));

depth = channel_positions(primary_channels+1 , 2); % want y coordinte
% does not match exactly...

%% Filter quality spikes
% by manually labeled
% [cids, cgs] = readClusterGroupsCSV(fullfile(phy_folder, 'cluster_group.tsv'));
% % or
% cids = T.cluster_id;
% cgs = T.group; % need to map to id?

% by quality metrics
% Thresholds to use: fr>0.01 && snr>4.0 && isi_violations_ratio<1
goodid = T.snr>4.0 & T.fr>0.01 & T.isi_violations_ratio<1;

goodid = (T.snr > 2.0 & T.fr > 0.05 & T.nn_hit_rate > 0.5 & ...
    T.isi_violations_ratio<1 & T.amplitude_cutoff<0.1 & T.presence_ratio>0.9); % Allen + some extra from isabel

%% Align spikes to teensy with TTLs
% honestly, this might be better using my python script

% 0) sync with existing TTL data
load(TTL_folder);
%spike_times_teensy = polyval(NIDQtoTeensy, polyval(NPtoNIDQ, spike_times));

spike_times_NIDQ = NPtoNIDQ(1)*double(spike_times) + NPtoNIDQ(2);
spike_times_teensy = NIDQtoTeensy(1)*spike_times_NIDQ + NIDQtoTeensy(2);

%% Test if syncing works here with beh data?
try
choice_teensy = [info.choice_time{:}];
catch
    tmp = cellfun(@(v) v(1), info.choice_time, 'un', 0);
    choice_teensy = [tmp{:}];
end

choice_nidq = nidq_choice1;
choice_nidq_toTeensy = NIDQtoTeensy(1) * double(choice_nidq) + NIDQtoTeensy(2);

%% Reformat spikes by good unit
%** this is indexing by 1 instead of 0!
spike_times_byunit = cell(sum(goodid),1);
goodunit = find(goodid);
for i = 1:length(goodunit)
    spike_times_byunit{i} = spike_times_teensy(spike_clusters==(goodunit(i)-1));    % fixed indexing!
end
try
    unit_depth = T.depth(goodunit); % from phy
catch
    unit_depth = depth(goodunit); % from kilosort
end


%% Redo making tons of plots, but with the new plotting code


unit = 23; % filter by depth as well?
dosave = 0;

savepath = 'D:\ANALYSIS\Kevin\Plots\KM47\250815\units_by_choice\';
if ~exist(savepath); mkdir(savepath); end


% for unit = 1:length(goodunit)

spiketimes = milliseconds(spike_times_byunit{unit});
% trialStarts = milliseconds([info.choice_time{:}]);
 [choice, reward, choice_time, ~, ~, ~, LTProb, ~] ...
    = extract_session_params(info, 1, 1);
trialStarts = milliseconds(choice_time);
% spiketimes_all = []; triallabels_all = [];
% for trial = 1:length(trialStarts)
%     spiketimes_all = [spiketimes_all, spiketimes' - trialStarts(trial)];
%     triallabels_all = [triallabels_all, trial * ones(1,length(spiketimes))];
% end

% g = figure;
% s = spikeRasterPlot(spiketimes_all,triallabels_all);
% s.XLimits = seconds([-5, 10]);

spiketimes_all = {}; % Mx1 cell array of spike times
for trial = 1:length(trialStarts)
    suse = spiketimes' - trialStarts(trial);
    suse(abs(suse)>seconds(60)) = []; % lets not save every spike every time
    spiketimes_all{trial} = seconds(  suse  );
end

g = figure;
plotSpikeRaster(spiketimes_all,'PlotType','vertline','XLimForCell',[-4 6]);
xlabel('Time relative to choice (s)')

if dosave
    u = goodunit(unit);
    spath = fullfile(savepath, [num2str(u) '--' 'depth' num2str(unit_depth(unit)) '--alltrials.png']);
    saveas(g, spath);
    close(g);
end





g=figure; hold on;
bin_centers = bin_edges(1:end-1)+seconds(bin_size/2);
plot(bin_centers, R_firing_rates)
plot(bin_centers, U_firing_rates)
xlabel('Time (s)') ; ylabel('Firing Rate (hz)');
legend('Reward','Unreward')



% by choice and reward...
% choice = [info.choice{:}]>=2;
% reward = [info.reward{:}];
choice = choice>=2;
g = figure;

% right reward
subplot(2,2,1);
trialID = find(choice==0 & reward==1); % Right reward
plotSpikeRaster(spiketimes_all(trialID),'PlotType','vertline','XLimForCell',[-4 6]);
title('Right reward')
% generate psth
bin_size = 0.1;
bin_edges = (-5:bin_size:10);
spike_counts = histcounts(horzcat(spiketimes_all{trialID}), bin_edges);
num_trials = length(trialID);
RR_firing_rates = spike_counts / (num_trials * bin_size); 



% right unreward
subplot(2,2,2);
trialID = find(choice==0 & reward==0); % Right unreward
plotSpikeRaster(spiketimes_all(trialID),'PlotType','vertline','XLimForCell',[-4 6]);
title('Right unreward')
% generate psth
bin_size = 0.1;
bin_edges = (-5:bin_size:10);
spike_counts = histcounts(horzcat(spiketimes_all{trialID}), bin_edges);
num_trials = length(trialID);
RU_firing_rates = spike_counts / (num_trials * bin_size); 

% left reward
subplot(2,2,3);
trialID = find(choice==1 & reward==1); % Left reward
plotSpikeRaster(spiketimes_all(trialID),'PlotType','vertline','XLimForCell',[-4 6]);
title('Left reward')
% generate psth
bin_size = 0.1;
bin_edges = (-5:bin_size:10);
spike_counts = histcounts(horzcat(spiketimes_all{trialID}), bin_edges);
num_trials = length(trialID);
LR_firing_rates = spike_counts / (num_trials * bin_size); 

% left reward
subplot(2,2,4);
trialID = find(choice==1 & reward==0); % Left reward
plotSpikeRaster(spiketimes_all(trialID),'PlotType','vertline','XLimForCell',[-4 6]);
title('Left unreward')
% generate psth
bin_size = 0.1;
bin_edges = (-5:bin_size:10);
spike_counts = histcounts(horzcat(spiketimes_all{trialID}), bin_edges);
num_trials = length(trialID);
LU_firing_rates = spike_counts / (num_trials * bin_size); 


if dosave
    spath = fullfile(savepath, [num2str(u) '--' 'depth' num2str(unit_depth(unit)) '-bychoice.png']);
    saveas(g, spath);
    close(g);
end

%end

% Plot average firing rates of every unit for each thing



g=figure; hold on;
bin_centers = bin_edges(1:end-1)+seconds(bin_size/2);
plot(bin_centers, RR_firing_rates)
plot(bin_centers, RU_firing_rates)
plot(bin_centers, LR_firing_rates)
plot(bin_centers, LU_firing_rates)
xlabel('Time (s)') ; ylabel('Firing Rate (hz)');
legend('Right reward','Right unreward','Left reward','Left unreward')

if dosave
    spath = fullfile(savepath, [num2str(u) '--' 'depth' num2str(unit_depth(unit)) '--firingrates.png']);
    saveas(g, spath);
    close(g);
end

% end

% neuron 66 is selective?

%% * old * 
% generate psth
bin_size = 0.1;
bin_edges = (-5:bin_size:10);
spike_counts = histcounts(horzcat(spiketimes_all{trialID}), bin_edges);
num_trials = length(trialID);
RR_firing_rates = spike_counts / (num_trials * bin_size); 



% pull and align all spike times
spiketimes = milliseconds(spike_times_byunit{unit});
trialStarts = milliseconds([info.choice_time{:}]);
spiketimes_all = []; triallabels_all = [];
for trial = 1:length(trialStarts)
    spiketimes_all = [spiketimes_all, spiketimes' - trialStarts(trial)];
    triallabels_all = [triallabels_all, trial * ones(1,length(spiketimes))];
end

% Generate firing rates for each condition

% subselect trials
trialID = find(choice==0 & reward==1); % Right reward
trial_subset = ismember(triallabels_all, trialID);
% generate psth
bin_size = 0.1;
bin_edges = seconds(-5:bin_size:10);
spike_counts = histcounts(spiketimes_all(trial_subset), bin_edges);
num_trials = length(trialID);
RR_firing_rates = spike_counts / (num_trials * bin_size); 

trialID = find(choice==0 & reward==0); % Right unreward
trial_subset = ismember(triallabels_all, trialID);
% generate psth
bin_size = 0.1;
bin_edges = seconds(-5:bin_size:10);
spike_counts = histcounts(spiketimes_all(trial_subset), bin_edges);
num_trials = length(trialID);
RU_firing_rates = spike_counts / (num_trials * bin_size); 

trialID = find(choice==1 & reward==1); % Left reward
trial_subset = ismember(triallabels_all, trialID);
% generate psth
bin_size = 0.1;
bin_edges = seconds(-5:bin_size:10);
spike_counts = histcounts(spiketimes_all(trial_subset), bin_edges);
num_trials = length(trialID);
LR_firing_rates = spike_counts / (num_trials * bin_size); 

trialID = find(choice==1 & reward==0); % Left unreward
trial_subset = ismember(triallabels_all, trialID);
% generate psth
bin_size = 0.1;
bin_edges = seconds(-5:bin_size:10);
spike_counts = histcounts(spiketimes_all(trial_subset), bin_edges);
num_trials = length(trialID);
LU_firing_rates = spike_counts / (num_trials * bin_size); 

%% test loading templates?
templates = readNPY(fullfile(phy_folder, 'templates.npy'));
templates_ind = readNPY(fullfile(phy_folder, 'template_ind.npy'));

tmp = squeeze(templates(200,:,:));
figure;plot(tmp)


%% Test regression model!

% want to predict spike count from : choice, reward, and value estimates of
% the current trial and one trial back
% - using a sliding window
% - as in previous papers

% todo:
% 1) Extract trial information
% 2) Extract value information from model fit
% 3) 

