%% Very similar test script, but for social

% IDK how good the session even is. I should make a plot...

% set paths
cd('D:\ANALYSIS\Kevin')
addpath(genpath('D:\ANALYSIS\Kevin\Utilities'))

% set folders
phy_folder = 'D:\DATA\Kevin\KM33-34\241231_g0\Spike_Sorting\phy';
TTL_folder = 'D:\DATA\Kevin\KM33-34\241231_g0\TTLs\TTLs.mat';
beh_folder = 'D:\DATA\Kevin\KM33-34\241231*.dat';
% mouse_name = 'KM33-34';
% from cluster
% phy_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM33\241219_g0\Spike_Sorting\phy';
% TTL_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM33\241219_g0\TTLs\TTLs.mat';
% beh_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM33\241219*.dat';
%  mouse_name = 'KM33-34';


phy_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM33-34\250111_g0\Spike_Sorting\phy';
TTL_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM33-34\250111_g0\TTLs\TTLs.mat';
beh_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM33-34\250111*.dat';
mouse_name = 'KM33-34';

% % 
phy_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM33-34\250128_g0\Spike_Sorting\phy';
TTL_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM33-34\250128_g0\TTLs\TTLs.mat';
beh_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM33-34\250128*.dat';
mouse_name = 'KM33-34';


phy_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM37-38\250213_g0\Spike_Sorting\phy';
TTL_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM37-38\250213_g0\TTLs\TTLs.mat';
beh_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM37-38\250213*.dat';
mouse_name = 'KM37-38';


phy_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM33-34\250208_g0\Spike_Sorting\phy';
TTL_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM33-34\250208_g0\TTLs\TTLs.mat';
beh_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM33-34\250208*.dat';
mouse_name = 'KM33-34';

phy_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM37-38\250224_g0\Spike_Sorting\phy';
TTL_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM37-38\250224_g0\TTLs\TTLs.mat';
beh_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM37-38\250224*.dat';
mouse_name = 'KM37-38';

phy_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM35-36\250225_g0\Spike_Sorting\phy';
TTL_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM35-36\250225_g0\TTLs\TTLs.mat';
beh_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM35-36\250225*.dat';
mouse_name = 'KM35-36';


phy_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM33-34\250218_g0\Spike_Sorting\phy';
TTL_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM33-34\250218_g0\TTLs\TTLs.mat';
beh_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM33-34\250218*.dat';
mouse_name = 'KM33-34';

phy_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM37-38\250228_g0\Spike_Sorting\phy';
TTL_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM37-38\250228_g0\TTLs\TTLs.mat';
beh_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM37-38\250228*.dat';
mouse_name = 'KM37-38';

phy_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM37-38\250211_g0\Spike_Sorting\phy';
TTL_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM37-38\250211_g0\TTLs\TTLs.mat';
beh_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM37-38\250211*.dat';
mouse_name = 'KM37-38';

phy_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM37-38\250314_g0\Spike_Sorting\phy';
TTL_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM37-38\250314_g0\TTLs\TTLs.mat';
beh_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM37-38\250314*.dat';
mouse_name = 'KM37-38';

phy_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM33-34\250318_g0\Spike_Sorting\phy';
TTL_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM33-34\250318_g0\TTLs\TTLs.mat';
beh_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM33-34\250318*.dat';
mouse_name = 'KM33-34';


phy_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM41-42\250420_g0\Spike_Sorting\phy';
TTL_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM41-42\250420_g0\TTLs\TTLs.mat';
beh_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM41-42\250420*.dat';
mouse_name = 'KM41-42';

% test the concatenated recordings?
%**TODO***
phy_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM33-34\250416_g0\Spike_Sorting_concatenated\phy';
TTL_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM33-34\250416_g0\TTLs\TTLs.mat';
beh_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM33-34\250416*.dat';
mouse_name = 'KM33-34';

% this is split in two. I did the syncing for, but just have to write some
% script to apply separately?
phy_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM41-42\250427_g0\Spike_Sorting_concatenated\phy';
TTL_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM41-42\250427_g0\TTLs\TTLs.mat';
TTL_folder_2 = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM41-42\250427-2_g0\TTLs\TTLs.mat';
beh_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM41-42\250427*.dat';
mouse_name = 'KM41-42';

phy_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM41-42\250429_g0\Spike_Sorting\phy';
TTL_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM41-42\250429_g0\TTLs\TTLs.mat';
beh_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM41-42\250429*.dat';
mouse_name = 'KM41-42';

% This is a social neutral session
phy_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM33-34\250408_g0\Spike_Sorting\phy';
TTL_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM33-34\250408_g0\TTLs\TTLs.mat';
beh_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM33-34\250408*.dat';
mouse_name = 'KM33-34';

% phy_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM37-38\250417_g0\Spike_Sorting\phy';
% TTL_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM37-38\250417_g0\TTLs\TTLs.mat';
% beh_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM37-38\250417*.dat';
% mouse_name = 'KM37-38';


phy_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM35-36\250414_g0\Spike_Sorting\phy';
TTL_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM35-36\250414_g0\TTLs\TTLs.mat';
beh_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM35-36\250414*.dat';
mouse_name = 'KM35-36';

% VCA1
phy_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM41-42\250728_g0\Spike_Sorting\phy';
TTL_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM41-42\250728_g0\TTLs\TTLs.mat';
beh_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM41-42\250728*.dat';
mouse_name = 'KM41-42';

% supposedly bla and stuff but i totally missed and ended up in a fiber
% tract
phy_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM43-44\250731_g0\Spike_Sorting\phy';
TTL_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM43-44\250731_g0\TTLs\TTLs.mat';
beh_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM43-44\250731*.dat';
mouse_name = 'KM43-44';

% also supposedly bla and ca2 and ppc but idk
phy_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM45-46\250802_g0\Spike_Sorting\phy';
TTL_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM45-46\250802_g0\TTLs\TTLs.mat';
beh_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM45-46\250802*.dat';
mouse_name = 'KM45-46';

% new vca1, but only have like 20 neurons idk why
phy_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM47-48\250809_g0\Spike_Sorting\phy';
TTL_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM47-48\250809_g0\TTLs\TTLs.mat';
beh_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM47-48\250809*.dat';
mouse_name = 'KM47-48';

% KM45 in BLA only
phy_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM45-46\250904_g0\Spike_Sorting\phy';
TTL_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM45-46\250904_g0\TTLs\TTLs.mat';
beh_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM45-46\250904*.dat';
mouse_name = 'KM45-46';
%%
% set paths
clc; clear; close all;
addpath(genpath('Z:\HarveyLab\Tier1\Kevin\Analysis\20250718_backup_Cindys_PC\Utilities'))

% 
% phy_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM33-34\250222_g0\Spike_Sorting\phy';
% TTL_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM33-34\250222_g0\TTLs\TTLs.mat';
% beh_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM33-34\250222*.dat';
% mouse_name = 'KM33-34';


phy_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM41-42\250805_g0\Spike_Sorting\phy';
TTL_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM41-42\250805_g0\TTLs\TTLs.mat';
beh_folder = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM41-42\250805*.dat';
mouse_name = 'KM41-42'; %**0 1 swap
%% Load data

% load unit info
spike_times = readNPY(fullfile(phy_folder, 'spike_times.npy'));
spike_clusters = readNPY(fullfile(phy_folder, 'spike_clusters.npy'));

% load auxiliary info <- this only works after running through phy
try
   
    T = readtable(fullfile(phy_folder, 'cluster_info.tsv'),...
        'FileType','text','Delimiter','\t');
    disp('loading quality metrics from phy')
% alternatively, all i need is depht, cluster_id, group, snr, fr,
% isi_violation
% 
% T_isiviolation = readtable(fullfile(phy_folder, 'cluster_isi_violations_ratio.tsv'),...
%      'FileType','text','Delimiter','\t');
% T_snr = readtable(fullfile(phy_folder, 'cluster_snr.tsv'),...
%      'FileType','text','Delimiter','\t');
% alternatively alternatively read from how saved
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


% load behavior data
beh_path = dir(beh_folder);
info = createBEHstruct_social(mouse_name, beh_path);

%%% Filter spikes by quality metrics
% Thresholds to use: fr>0.01 && snr>4.0 && isi_violations_ratio<1
goodid = T.snr>4.0 & T.fr>0.01 & T.isi_violations_ratio<1; % isabel
%goodid = T.snr>1.0 & T.fr>0.5 & T.isi_violations_ratio<1 & T.presence_ratio>0.5; % lily

goodid = T.snr>1.0 & T.fr>0.01 & T.isi_violations_ratio<1 & T.presence_ratio>0.5; % allen?

goodid = (T.snr > 2.0 & T.fr > 0.05 & T.nn_hit_rate > 0.5 & ...
    T.isi_violations_ratio<1 & T.amplitude_cutoff<0.1 & T.presence_ratio>0.9); % Allen + some extra from isabel

%goodid = (T.snr>1.0 & T.isi_violations_ratio<1 & T.amplitude_cutoff<0.1 & T.presence_ratio>0.9);



%%% Align spikes to teensy with TTLs
%*** this is not working for the social mice...
% honestly, this might be better using my python script

% 0) sync with existing TTL data
load(TTL_folder);
%spike_times_teensy = polyval(NIDQtoTeensy, polyval(NPtoNIDQ, spike_times));

spike_times_NIDQ = NPtoNIDQ(1)*double(spike_times) + NPtoNIDQ(2);
spike_times_teensy = NIDQtoTeensy(1)*spike_times_NIDQ + NIDQtoTeensy(2);

%%% Align spikes with TTLs for split files
% to test, need to see if offset subtractino works ok
if 0 % only run when concatenating
% piecewise
load(TTL_folder);
nCatThresh = nCat;
spike_times_NIDQ = NPtoNIDQ(1) * double(spike_times(spike_times<nCatThresh)) + NPtoNIDQ(2);
spike_times_teensy_1 = NIDQtoTeensy(1)*spike_times_NIDQ + NIDQtoTeensy(2);


load(TTL_folder_2);
% there needs to be some offset thing here... not sure if ive done it right
%** no I havent.
spike_times_NIDQ = NPtoNIDQ(1) * double(spike_times(spike_times>nCatThresh) - nCatThresh) + NPtoNIDQ(2) ;
spike_times_teensy_2 = NIDQtoTeensy(1)*spike_times_NIDQ + NIDQtoTeensy(2);

spike_times_teensy = [spike_times_teensy_1; spike_times_teensy_2];

%** use choice times to sync?
choice_time_all = vertcat(info.choice_time{:});
choice_time_nidq = sort([nidq_choice1, nidq_choice2]);
choice_time_nidq_to_teensy = NIDQtoTeensy(1) * choice_time_nidq;
% find match indexes?

% compute *new* sync line

% retry the syncs?
end

%%% Test if syncing works here with beh data?
% *** this wont work for sessions before 2025/01/07 since I'm not writing
% the heartbeat light for every choice (dumb mistake on my end)

%choice_teensy = [info.choice_time{:}];
% updated above for social with multipel choice times
mouseID = info.mouseID;
choice_time = info.choice_time; choice_time_inferred = info.choice_time_inferred;
choice_time_1 = []; choice_time_1_inferred = [];
choice_time_2 = []; choice_time_2_inferred = [];
for trial = 1:length(choice_time)
    for i = 1:length(choice_time{trial})
        if mouseID{trial}(i) == 0
            choice_time_1(end+1) = choice_time{trial}(i);
            choice_time_1_inferred(end+1) = choice_time_inferred{trial}(i);
        elseif mouseID{trial}(i) == 1
            choice_time_2(end+1) = choice_time{trial}(i);
            choice_time_2_inferred(end+1) = choice_time_inferred{trial}(i);
        end
    end
end

choice_nidq = nidq_choice1;
choice_nidq1_toTeensy = NIDQtoTeensy(1) * double(choice_nidq) + NIDQtoTeensy(2);
choice_nidq = nidq_choice2;
choice_nidq2_toTeensy = NIDQtoTeensy(1) * double(choice_nidq) + NIDQtoTeensy(2);

% these should match up within 1 ms
figure; plot(choice_nidq1_toTeensy - choice_time_1);
hold on; plot(choice_nidq2_toTeensy - choice_time_2);

%%% Preprocess for plotting, not just checking syncing?

% load data
choice = info.choice; choice_time = info.choice_time;
reward = info.reward;
choice_time_inferred = info.choice_time_inferred;
badid = cellfun(@isempty, choice) | cellfun(@isempty, reward) | ...
    cellfun(@length, choice) ~= cellfun(@length, reward);
choice(badid) = []; reward(badid) = [];
mouseID = info.mouseID; mouseID(badid) = [];
choice_time(badid) = []; choice_time_inferred(badid) = [];


% format data
choice1 = zeros(1,length(choice));  % left choice is 1, right is -1!
reward1 = zeros(1,length(choice)); 
choice2 = zeros(1,length(choice)); 
reward2 = zeros(1,length(choice)); 
choice_time_1 = zeros(1,length(choice));
choice_time_2 = zeros(1,length(choice));
choice_time_1_inferred = zeros(1,length(choice));
choice_time_2_inferred = zeros(1,length(choice));
outcomevals = [-1, 1]; % -1 = right, 1 = left
for trial = 1:length(choice)
    for i = 1:length(choice{trial})
        if mouseID{trial}(i) == 1 %0
            % set so assings 1 or -1, and unsassigned is auto 0?
            choice1(trial) = outcomevals( (double(choice{trial}(i)) > 1)+1 );
            reward1(trial) = outcomevals( reward{trial}(i)+1 );
            choice_time_1(trial) = choice_time{trial}(i);
            choice_time_1_inferred(trial) = choice_time_inferred{trial}(i);
        elseif mouseID{trial}(i) == 0 %1
            choice2(trial) = outcomevals( (double(choice{trial}(i)) > 1)+1 );
            reward2(trial) = outcomevals( reward{trial}(i)+1 );
            choice_time_2(trial) = choice_time{trial}(i);
            choice_time_2_inferred(trial) = choice_time_inferred{trial}(i);
        end
    end
end
% filter trials with one choice
badid = choice1==0 | choice2==0;
choice1(badid) = []; choice2(badid) = [];
reward1(badid) = []; reward2(badid) = [];
choice_time_1(badid) = []; choice_time_2(badid) = [];
choice_time_1_inferred(badid) = []; choice_time_2_inferred(badid) = [];

%% pull above, but also q values?
% [Qmulti,~,acc_multi, Qother] = extract_value_information(info, 'multi-agent');
% % this is a terrible fit

model_use = 'non-social';
[Q1,~, ~] = extract_value_information(info, model_use);
[Q2, ~, acc] = extract_value_information(info, model_use, 2); 

deltaQ = Q1(1,:) - Q1(2,:);
deltaQ2 = Q2(1,:) - Q2(2,:);
figure; plot(deltaQ, deltaQ2, 'o');

value_bins = -1:0.2:1;
value_centers = (value_bins(1:end-1) + value_bins(2:end))/2;

[counts1, ~, idx] = histcounts(deltaQ, value_bins);
[counts2, ~, idx2] = histcounts(deltaQ2, value_bins);
counts1
counts2

%% Reformat spikes, and filter by good units

spike_times_byunit = cell(sum(goodid),1);
goodunit = find(goodid);
for i = 1:length(goodunit)
    spike_times_byunit{i} = spike_times_teensy(spike_clusters==(goodunit(i)-1)); % fixed indexing!    
end

% these two dont match and I cannot figure out why
try
    unit_depth = T.depth(goodunit); % from phy
catch
    unit_depth = depth(goodunit); % from kilosort
end

%%% Filter for crashed recording
crash_time = max(spike_times_teensy);
badid_1 = find(choice_time_1 > crash_time);
badid_2 = find(choice_time_2 > crash_time);
crash_trial = min([badid_1, badid_2]);

badid = crash_trial : length(choice1);
choice1(badid) = []; choice2(badid) = [];
reward1(badid) = []; reward2(badid) = [];
choice_time_1(badid) = []; choice_time_2(badid) = [];
choice_time_1_inferred(badid) = []; choice_time_2_inferred(badid) = [];

%%% I might redo all the plotting code
% - since the code I downloaded DOES NOT let me modify *anything*
% - i had a great script I was using in grad school idk what happened.

unit = 50; % and 146 might be choice?

plotOtherMouse = 1;
plotAccoutrements = 1;
plotValue = 0;

spiketimes = milliseconds(spike_times_byunit{unit});
trialStarts = milliseconds(choice_time_1_inferred); % <- the first lick

spiketimes_all = {}; % Mx1 cell array of spike times
for trial = 1:length(trialStarts)
    suse = spiketimes' - trialStarts(trial);
    suse(abs(suse)>seconds(60)) = []; % lets not save every spike every time
    spiketimes_all{trial} = seconds(  suse  );
end

g = figure;
plotSpikeRaster(spiketimes_all,'PlotType','vertline','XLimForCell',[-4 10]);

% get other mouse?
if plotOtherMouse
    trialStartOther = milliseconds(choice_time_2_inferred); 
    trialStartOther = seconds( trialStartOther - trialStarts );
    hold on;
    plot( trialStartOther, (1:length(trialStarts)), 'ro','MarkerFaceColor','r','MarkerSize',3);
end
if plotAccoutrements % this assumes everythign is aligned to first lick!
    trialReward = milliseconds( choice_time_1 );
    trialReward = seconds(trialReward - trialStarts);
    plot(trialReward, (1:length(trialStarts)), 'bx','MarkerFaceColor', 'b','MarkerSize',5);

    trialOtherReward = milliseconds( choice_time_2 );
    trialOtherReward = seconds(trialOtherReward - trialStarts);
    plot(trialOtherReward,  (1:length(trialStarts)), 'rx','MarkerFaceColor','r','MarkerSize',5);
    
    % plot choice info?
    tmp = xlim;
    plot(tmp(1)*ones(1,sum(choice1==-1 & reward1==1)), find(choice1==-1 & reward1==1)  , 'b>', 'MarkerFaceColor','b', 'MarkerSize',5) % right
    plot(tmp(1)*ones(1,sum(choice1==-1 & reward1==-1)), find(choice1==-1 & reward1==-1)  , 'b>', 'MarkerSize',5) % right
    plot(0.25+(tmp(1)*ones(1,sum(choice1==1 & reward1==1))), find(choice1==1 & reward1==1)  , 'b<', 'MarkerFaceColor','b', 'MarkerSize',5) % right
    plot(0.25+(tmp(1)*ones(1,sum(choice1==1 & reward1==-1))), find(choice1==1 & reward1==-1)  , 'b<', 'MarkerSize',5) % right

    plot(-1+(tmp(1)*ones(1,sum(choice2==-1 & reward2==1))), find(choice2==-1 & reward2==1)  , 'r>', 'MarkerFaceColor','r', 'MarkerSize',5) % right
    plot(-1+(tmp(1)*ones(1,sum(choice2==-1 & reward2==-1))), find(choice2==-1 & reward2==-1)  , 'r>', 'MarkerSize',5) % right
    plot(-0.75+(tmp(1)*ones(1,sum(choice2==1 & reward2==1))), find(choice2==1 & reward2==1) , 'r<', 'MarkerFaceColor','r', 'MarkerSize',5) % right
    plot(-0.75+(tmp(1)*ones(1,sum(choice2==1 & reward2==-1))), find(choice2==1 & reward2==-1)  , 'r<', 'MarkerSize',5) % right
end

xlabel('Time of choice (port lick)')

%% new plot everything
% lets use the new raster plot
% and also plot more controlled
% e.g., plot right reward, given 4 other conditions
%       and left reward, given 4 other conditions

% maybe also plot when the othe rmouse licked?

% the session data MUST be processed in a cell above
% this needs choice_time_1, 
%% align to choice all
dosave = 0;
% unit=120;
disp('make sure this is set up for the right mouse in social!!!')

savepath = 'D:\ANALYSIS\Kevin\Plots\KM45-46\250904\socialplus\';
if ~exist(savepath); mkdir(savepath); end

for unit = 1:length(goodunit)

%
unit
spiketimes = milliseconds(spike_times_byunit{unit});
trialStarts = milliseconds(choice_time_1); % <- when the choice detectino occurs
trialStarts = milliseconds(choice_time_1_inferred); % <- the first lick

trialStartOther = milliseconds(choice_time_2_inferred); 
trialStartOther = seconds( trialStartOther - trialStarts );

% extract spike times
spiketimes_all_struct.(strcat("unit",string(unit)))= {}; % Mx1 cell array of spike times
for trial = 1:length(trialStarts)
    suse = spiketimes' - trialStarts(trial);
    suse(abs(suse)>seconds(60)) = []; % lets not save every spike every time
    spiketimes_all_struct.(strcat("unit",string(unit))){trial} = seconds(  suse  );
end
end
% plot all spikes, but were actually just going to do this another way
% g = figure;
% plotSpikeRaster(spiketimes_all,'PlotType','vertline','XLimForCell',[-4 6]);
% xlabel('Time relative to choice (s)')
%
% left reward self

%% to firing rate (binning, also bin choice time, reward times)

bin_size = 0.1;
bin_edges = (-4:bin_size:6);
for unit = 1:length(goodunit)
FR_all.(strcat("unit",string(unit))) = [];
end
binned_trialReward=[];

binned_trialStartOther=[];

binned_trialOtherReward=[];
for trial = 1:length(trialStarts)
    trial
for unit = 1:length(goodunit)
%


FR_all.(strcat("unit",string(unit)))(:,trial) = histcounts(spiketimes_all_struct.(strcat("unit",string(unit))){trial}, bin_edges);
end
binned_trialReward(:,trial)=histcounts(trialReward(trial), bin_edges);

binned_trialStartOther(:,trial)=histcounts(trialStartOther(trial), bin_edges);

binned_trialOtherReward(:,trial)=histcounts(trialOtherReward(trial), bin_edges);
end
%% population dynamics in single trials
%reformat single unit concatenated single trial
units = fieldnames(FR_all);         % {'unit1','unit2',...}
nUnits = numel(units);
nTrials = size(FR_all.(units{1}),2);

FR_concat=[];

for tr = 1:nTrials
    % collect all units for this trial
    trialData = [];
    for u = 1:nUnits
        trialData(:,u) = FR_all.(units{u})(:,tr);
    end
    FR_concat = [FR_concat;zscore(trialData,0,1)];  %z-scored
    % each field = (timeBins × nUnits)
end
%% PCA/Umap
dim_red="PCA"; %"Umap"
if dim_red=="PCA"
[coeff,score,~,~,explained,~] = pca(FR_concat);
elseif dim_red=="Umap"
n_neighbors=50;
addpath(genpath('umap'));
[score, umap] = run_umap(FR_concat, 'n_components', 2,'n_neighbors',n_neighbors, 'verbose', 'none','randomize',true);
end


%% four conditions
timeWithinConcat = repmat(bin_edges(2:end), 1, nTrials);
trialIDMat = repmat([1:nTrials],100,1);
trialIDConcat = reshape(trialIDMat, [],1); %trial ID in concatenated series

binned_trialStartOtherConcat=reshape(binned_trialStartOther, [],1);
binned_trialRewardConcat=reshape(binned_trialReward, [],1);
span=30
dimension=2
condition_c1=[-1 1];
condition_c2=[-1 1];
condition_r1=[-1 1];
condition_r2=[-1 1];
for pc= 1:dimension %PC
    score_sm(:,pc)=smooth(score(:,pc),span,'sgolay',3);
end

range=[-1.5 2];
% figure;
% bar(explained(1:10))
% hold on
% plot(cumsum(explained(1:10)),'-o')
i=1;
color_seq = orderedcolors("gem");
figure;

labels=[]
subplot(10,1,[1 6])
 scatter(score_sm(:,1),score_sm(:,2),30,timeWithinConcat,"MarkerEdgeAlpha",0.2)
 hold on
 for c1=condition_c1
     for c2=condition_c2
         % for r1=condition_r1 %condition_r1
         %    for  r2=condition_r2%condition_r2
         trial_ID=find(choice1==c1 & choice2==c2);
         % trial_ID=find(reward1==r1 & reward2==r2);
         [tf, condition_time] = ismember(trialIDConcat, trial_ID);
         subplot(10,1,[1 6])
         % scatter(score_sm(find(timepointsConcat),1),score_sm(find(timepointsConcat),2),[],timeWithinConcat(find(timepointsConcat)),"MarkerEdgeAlpha",0.3)
         hold on
         % plot(score_sm(find(condition_time),1),score_sm(find(condition_time),2), "LineWidth",1,"Color",[color_seq(i,:) 0.07])
         pc1=score_sm(find(condition_time),1);
         pc2=score_sm(find(condition_time),2);
          % scatter(score_sm(find(condition_time' & timeWithinConcat==0),1),score_sm(find(condition_time' & timeWithinConcat==0),2),30,"filled", ...
         %     "MarkerFaceColor",color_seq(i,:),"MarkerEdgeColor","w","Marker","^","MarkerFaceAlpha",0.7)
         scatter(score_sm(find(condition_time' & binned_trialStartOtherConcat'==1),1),score_sm(find(condition_time' & binned_trialStartOtherConcat'==1),2),30,"filled", ...
             "MarkerFaceColor",color_seq(i,:),"MarkerEdgeColor","w","Marker","o","MarkerFaceAlpha",0.7)
       
         m=plot(mean(reshape(pc1,100,[]),2),mean(reshape(pc2,100,[]),2), "LineWidth",4,"Color",color_seq(i,:),'DisplayName',strcat("c1:",string(c1), " c2:",string(c2)))
         scatter(mean(score_sm(find(condition_time' & timeWithinConcat==0),1)),mean(score_sm(find(condition_time' & timeWithinConcat==0),2)),80,"filled", ...
             "MarkerFaceColor",color_seq(i,:),"MarkerEdgeColor","k","Marker","^")
          % xlim(range)
         subplot(10,1,i+6)
         histogram(score_sm(find(condition_time' & binned_trialStartOtherConcat'==1),1),[range(1):0.2:range(2)],"FaceColor",color_seq(i,:))
         % xlim(range)
         i=i+1;

     end
end

%% two conditions
timeWithinConcat = repmat(bin_edges(2:end), 1, nTrials);
trialIDMat = repmat([1:nTrials],100,1);
trialIDConcat = reshape(trialIDMat, [],1); %trial ID in concatenated series

binned_trialStartOtherConcat=reshape(binned_trialStartOther, [],1);
binned_trialRewardConcat=reshape(binned_trialReward, [],1);
span=30
dimension=2
condition_c1=[-1 1];
condition_c2=[-1 1];
condition_r1=[-1 1];
condition_r2=[-1 1];

condition_label="self choice";
condition_def=choice2==condition_c2(1);


for pc= 1:dimension %PC
    score_sm(:,pc)=smooth(score(:,pc),span,'sgolay',3);
end

range=[-1.5 2];
% figure;
% bar(explained(1:10))
% hold on
% plot(cumsum(explained(1:10)),'-o')
i=1;
color_seq = orderedcolors("gem");
figure;

labels=[]
 % scatter(score_sm(:,1),score_sm(:,2),30,timeWithinConcat,"MarkerEdgeAlpha",0.2)
 hold on
 for c=[1 0]
     % for r1=condition_r1 %condition_r1
     %    for  r2=condition_r2%condition_r2
     trial_ID=find(condition_def==c);
     % trial_ID=find(reward1==r1 & reward2==r2);
     [tf, condition_time] = ismember(trialIDConcat, trial_ID);

     % scatter(score_sm(find(timepointsConcat),1),score_sm(find(timepointsConcat),2),[],timeWithinConcat(find(timepointsConcat)),"MarkerEdgeAlpha",0.3)
     hold on
     pc1=score_sm(find(condition_time),1);
     pc2=score_sm(find(condition_time),2);
     plot(pc1,pc2, "LineWidth",1,"Color",[color_seq(i,:) 0.1])

     % scatter(score_sm(find(condition_time' & timeWithinConcat==0),1),score_sm(find(condition_time' & timeWithinConcat==0),2),30,"filled", ...
     %     "MarkerFaceColor",color_seq(i,:),"MarkerEdgeColor","w","Marker","^","MarkerFaceAlpha",0.7)
     scatter(score_sm(find(condition_time' & binned_trialStartOtherConcat'==1),1),score_sm(find(condition_time' & binned_trialStartOtherConcat'==1),2),30,"filled", ...
         "MarkerFaceColor",color_seq(i,:),"MarkerEdgeColor","w","Marker","o","MarkerFaceAlpha",0.7)

     m=plot(mean(reshape(pc1,100,[]),2),mean(reshape(pc2,100,[]),2), "LineWidth",4,"Color",color_seq(i,:));
     scatter(mean(score_sm(find(condition_time' & timeWithinConcat==0),1)),mean(score_sm(find(condition_time' & timeWithinConcat==0),2)),80,"filled", ...
         "MarkerFaceColor",color_seq(i,:),"MarkerEdgeColor","k","Marker","^")
     % xlim(range)
     i=i+1;

 end


%% trial_average under conditions

condition_c1=[-1 1];
condition_c2=[-1 1];
condition_r1=[-1 1];
condition_r2=[-1 1];
FR_trial_average=[];
condition=[];
i=1;
for c1=condition_c1
    for c2=condition_c2
        for r1=condition_r1
            for r2=condition_r2
                condition=[condition [c1;c2;r1;r2]];
                trial_ID=find(choice1==c1 & choice2==c2 & reward1==r1 & reward2==r2);
                trial_count(i,:)= numel(trial_ID);
                for unit = 1:length(goodunit)
                    FR_trial_average.(strcat("unit",string(unit)))(:,i)=mean(FR_all.(strcat("unit",string(unit)))(:,trial_ID),2);

                end
                i=i+1;
            end
        end
    end
end
%% stackedplot
unit=15;
figure;
s = heatmap(FR_trial_average.(strcat("unit",string(unit)))',"GridVisible","off")
colormap("turbo")
h.XDisplayLabels = {string(bin_edges)};
% h.YDisplayLabels = {'A','B','C','D','E'};
%% single cell plotting
g=figure; 
count = 1;
c1 = -1; r1 = 1; % fixing on one mouse
r1 = 1; r2 = 1;% fixing on reward
unit=1;
psth_all = {}; label_all = {}; num_trials_all = [];
bin_size = 0.1;
bin_edges = (-5:bin_size:10);
for c1 = [-1,1] % change these as needed
    for c2 = [-1,1]
        % reward fixed
        %trialID = find(choice1==c1 & reward1==r1 & choice2==c2 & reward2==r2);
        % reward agnostic
        trialID = find(choice1==c1 & choice2==c2);
        subplot(2,2,count); hold on;
        plotSpikeRaster(spiketimes_all_struct.(strcat("unit",string(unit)))(trialID),'PlotType','vertline','XLimForCell',[4 6]);
        title(sprintf('c1: %d, c2: %d', c1,c2));
        label_all{count} = sprintf('c1: %d, c2: %d', c1,c2);
        % plot other mouse lick?
        hold on;
        plot( trialStartOther(trialID), (1:length(trialID)), 'ro','MarkerFaceColor','r','MarkerSize',3);
        % get psth data for later
        spike_counts = histcounts(horzcat(spiketimes_all_struct.(strcat("unit",string(unit))){trialID}), bin_edges);
        num_trials = length(trialID);
        psth_all{count} = spike_counts / (num_trials * bin_size); 
        num_trials_all(count) = num_trials;
        count = count+1;
    end
end



if dosave
    u = goodunit(unit);
    spath = fullfile(savepath, [num2str(u) '--' 'depth' num2str(unit_depth(unit)) '--raster.png']);
    saveas(g, spath);
    close(g);
end

% plot the other things
g=figure; hold on;
bin_centers = bin_edges(1:end-1)+seconds(bin_size/2);
axlabel=[];
for c = 1:length(psth_all)
    if num_trials_all(c) > 15
        plot(bin_centers, psth_all{c})
        axlabel(end+1) = c;
    end
end
xlabel('Time (s)') ; ylabel('Firing Rate (hz)');
legend(label_all(axlabel))

if dosave
    u = goodunit(unit);
    spath = fullfile(savepath, [num2str(u) '--' 'depth' num2str(unit_depth(unit)) '--psths.png']);
    saveas(g, spath);
    close(g);
end
% same/diff choice, reward/un-reward
g=figure; 
count = 1;
c1 = -1; r1 = 1; % fixing on one mouse
r1 = 1; r2 = 1;% fixing on reward

psth_all = {}; label_all = {}; num_trials_all = [];
bin_size = 0.1;
bin_edges = (-5:bin_size:10);
for c = [-1,1] % 1:same, -1:diff
    for r1 = [-1,1] % mouse1 rew/unrew
        % reward fixed
        %trialID = find(choice1==c1 & reward1==r1 & choice2==c2 & reward2==r2);
        % reward agnostic
        trialID = find(choice1.*choice2==c & reward1==r1);
        subplot(2,2,count); hold on;
        plotSpikeRaster(spiketimes_all_struct.(strcat("unit",string(unit)))(trialID),'PlotType','vertline','XLimForCell',[-4 6]);
        title(sprintf('choice: %d, r1: %d', c, r1));
        label_all{count} = sprintf('choice: %d, r1: %d', c,r1);
        % plot other mouse lick?
        hold on;
        plot( trialStartOther(trialID), (1:length(trialID)), 'ro','MarkerFaceColor','r','MarkerSize',3);
        % get psth data for later
        spike_counts = histcounts(horzcat(spiketimes_all_struct.(strcat("unit",string(unit))){trialID}), bin_edges);
        num_trials = length(trialID);
        psth_all{count} = spike_counts / (num_trials * bin_size); 
        num_trials_all(count) = num_trials;
        count = count+1;
    end
end


% plot the other things
g=figure; hold on;
bin_centers = bin_edges(1:end-1)+seconds(bin_size/2);
axlabel=[];
for p = 1:length(psth_all)
    if num_trials_all(p) > 15
        plot(bin_centers, psth_all{p})
        axlabel(end+1) = p;
    end
end
xlabel('Time (s)') ; ylabel('Firing Rate (hz)');
legend(label_all(axlabel))

%% Temp, figures for a grant

if 1
try
g=figure; 
count = 1;

psth_all = {}; num_trials_all = [];
bin_size = 0.1;
bin_edges = (-5:bin_size:10);
yticks_cat_save = [];

ncats = 0;
minN = 12;
cuse = linspecer(4);

ax1 = subplot(3,1, [1,2]);

for c1 = [-1,1] % change these as needed
    for c2 = [-1,1]
        % reward agnostic
        trialID = find(choice1==c1 & choice2==c2);
        % subsample
        trialID = trialID(randsample(1:length(trialID), minN));
        % add offset
        offset_cell = repmat({-60}, 1, ncats);
        spiketimes_plot = [offset_cell, spiketimes_all(trialID)];
        
        % plot
        LineFormat = struct();LineFormat.Color = cuse(count,:);

        plotSpikeRaster(spiketimes_plot,'PlotType','vertline','XLimForCell',[-4 6],'LineFormat',LineFormat);
        
        yticks_cat_save(count) = ncats + round(length(trialID)/2);

        % get psth data for later
        spike_counts = histcounts(horzcat(spiketimes_all{trialID}), bin_edges);
        num_trials = length(trialID);
        psth_all{count} = spike_counts / (num_trials * bin_size); 
        num_trials_all(count) = num_trials;

        % tabulate 
        count = count+1;
        ncats = ncats + length(trialID);
    end
end
yticks(yticks_cat_save)
yticklabels({'Both Left','Self left, other right','Self right, other left','Both right'})
xlabel('Time to lick (s)');


ax2 = subplot(3,1, 3); hold on;
%bin_centers = bin_edges(1:end-1)+seconds(bin_size/2);
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end))/2;
axlabel=[];
for c = 1:length(psth_all)
    if num_trials_all(c) > 10
        plot(bin_centers, psth_all{c}, 'Color', cuse(c,:));
        axlabel(end+1) = c;
    end
end
xlabel('Time (s)') ; ylabel('Firing Rate (hz)');

linkaxes([ax1, ax2], 'x');
if dosave
    u = goodunit(unit);
    spath = fullfile(savepath, [num2str(u) '--' 'depth' num2str(unit_depth(unit)) '--bychoice.png']);
    saveas(g, spath);
    close(g);
end


catch
    disp('skipping')
    close(g);
end

end

%% Temp value fig!
% - plot neuron sorted by binned value?
value_bins = -1:0.2:1;
value_centers = (value_bins(1:end-1) + value_bins(2:end))/2;
% 
deltaQ = Q1(1,:) - Q1(2,:);
[counts1, ~, idx] = histcounts(deltaQ, value_bins);

deltaQ2 = Q2(1,:) - Q2(2,:);
[counts2, ~, idx2] = histcounts(deltaQ2, value_bins);

% SAMPLE FROM THE ABOVE BINS
g=figure; 
count = 1;

psth_all = {}; num_trials_all = [];
bin_size = 0.1;
bin_edges = (-5:bin_size:10);
yticks_cat_save = [];

ncats = 0;
minN = 15; % num samps per bin
cuse = jet(length(value_centers));
%cuse = redblue(length(value_centers));
cuse = hot(round(length(value_centers)*1.5)); cuse = cuse(1:length(value_centers), :);

ax1 = subplot(3,2, [1,3]);

for vbin = 1:length(value_centers)
    if counts1(vbin) < minN; continue; end

    trialID = find(idx==vbin);
    trialID = trialID(randsample(1:length(trialID), minN));
    

    offset_cell = repmat({-60}, 1, ncats);
    spiketimes_plot = [offset_cell, spiketimes_all(trialID)];
        
    % plot
    LineFormat = struct();LineFormat.Color = cuse(count,:);

    plotSpikeRaster(spiketimes_plot,'PlotType','vertline','XLimForCell',[-4 6],'LineFormat',LineFormat);
        
    yticks_cat_save(count) = ncats + round(length(trialID)/2);

    % get psth data for later
    spike_counts = histcounts(horzcat(spiketimes_all{trialID}), bin_edges);
    num_trials = length(trialID);
    psth_all{count} = spike_counts / (num_trials * bin_size); 
    num_trials_all(count) = num_trials;

    % tabulate 
    count = count+1;
    ncats = ncats + length(trialID);
end
yticks(yticks_cat_save)
yticklabels(value_centers); ylabel('Relative value');
xlabel('Time to lick (s)');
title('Sorted by self value')


ax2 = subplot(3,2, 5); hold on;
%bin_centers = bin_edges(1:end-1)+seconds(bin_size/2);
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end))/2;
axlabel=[];
for c = 1:length(psth_all)
    if num_trials_all(c) > 4
        plot(bin_centers, psth_all{c}, 'Color', cuse(c,:));
        axlabel(end+1) = c;
    end
end
xlabel('Time (s)') ; ylabel('Firing Rate (hz)');

linkaxes([ax1, ax2], 'x');

% Other
count = 1;
psth_all = {}; num_trials_all = [];
yticks_cat_save = [];
ncats = 0;

ax3 = subplot(3,2, [2,4]);

for vbin = 1:length(value_centers)
    if counts2(vbin) < minN; continue; end
    trialID = find(idx2==vbin);
    trialID = trialID(randsample(1:length(trialID), minN));
    

    offset_cell = repmat({-60}, 1, ncats);
    spiketimes_plot = [offset_cell, spiketimes_all(trialID)];
        
    % plot
    LineFormat = struct();LineFormat.Color = cuse(count,:);

    plotSpikeRaster(spiketimes_plot,'PlotType','vertline','XLimForCell',[-4 6],'LineFormat',LineFormat);
        
    yticks_cat_save(count) = ncats + round(length(trialID)/2);

    % get psth data for later
    spike_counts = histcounts(horzcat(spiketimes_all{trialID}), bin_edges);
    num_trials = length(trialID);
    psth_all{count} = spike_counts / (num_trials * bin_size); 
    num_trials_all(count) = num_trials;

    % tabulate 
    count = count+1;
    ncats = ncats + length(trialID);
end
yticks(yticks_cat_save)
yticklabels(value_centers); ylabel('Relative value');
xlabel('Time to lick (s)');
title('Sorted by other value')


ax4 = subplot(3,2, 6); hold on;
%bin_centers = bin_edges(1:end-1)+seconds(bin_size/2);
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end))/2;
axlabel=[];
for c = 1:length(psth_all)
    if num_trials_all(c) > 4
        plot(bin_centers, psth_all{c}, 'Color', cuse(c,:));
        axlabel(end+1) = c;
    end
end
xlabel('Time (s)') ; ylabel('Firing Rate (hz)');

linkaxes([ax3, ax4], 'x');


if dosave
    u = goodunit(unit);
    spath = fullfile(savepath, [num2str(u) '--' 'depth' num2str(unit_depth(unit)) '--byvalue.png']);
    saveas(g, spath);
    close(g);
end


%%

% end

%% Same as above, but for not just rewarded trials?
savepath = 'D:\ANALYSIS\Kevin\Plots\KM33-34\250130\units_social\';
if ~exist(savepath); mkdir(savepath); end

for unit = 1:length(goodunit)

%%

spiketimes = milliseconds(spike_times_byunit{unit});
trialStarts = milliseconds(choice_time_1); % <- when the choice detectino occurs
trialStarts = milliseconds(choice_time_1_inferred); % <- the first lick

trialStartOther = milliseconds(choice_time_2_inferred); 
trialStartOther = seconds( trialStartOther - trialStarts );

% extract spike times
spiketimes_all = {}; % Mx1 cell array of spike times
for trial = 1:length(trialStarts)
    suse = spiketimes' - trialStarts(trial);
    suse(abs(suse)>seconds(60)) = []; % lets not save every spike every time
    spiketimes_all{trial} = seconds(  suse  );
end

% plot all spikes, but were actually just going to do this another way
% g = figure;
% plotSpikeRaster(spiketimes_all,'PlotType','vertline','XLimForCell',[-4 6]);
% xlabel('Time relative to choice (s)')

% left reward self
g=figure; 
count = 1;

psth_all = {}; label_all = {}; num_trials_all = [];
bin_size = 0.1;
bin_edges = (-5:bin_size:10);
for r1 = [-1,1] % change these as needed
    for r2 = [-1,1]
        for c1 = [-1,1]
            for c2 = [-1,1]
                trialID = find(choice1==c1 & reward1==r1 & choice2==c2 & reward2==r2);
                subplot(4,4,count); hold on;
                plotSpikeRaster(spiketimes_all(trialID),'PlotType','vertline','XLimForCell',[-4 6]);
                title(sprintf('c1: %d, c2: %d, r1: %d, r2: %d', c1,c2,r1,r2));
                label_all{count} = sprintf('c1: %d, c2: %d, r1: %d, r2: %d', c1,c2,r1,r2);
                % plot other mouse lick?
                hold on;
                plot( trialStartOther(trialID), (1:length(trialID)), 'ro','MarkerFaceColor','r','MarkerSize',3);
                % get psth data for later
                spike_counts = histcounts(horzcat(spiketimes_all{trialID}), bin_edges);
                num_trials = length(trialID);
                psth_all{count} = spike_counts / (num_trials * bin_size); 
                num_trials_all(count) = num_trials;
                count = count+1;
            end
        end
    end
end



if dosave
    u = goodunit(unit);
    spath = fullfile(savepath, [num2str(u) '--' 'depth' num2str(unit_depth(unit)) '--raster.png']);
    saveas(g, spath);
    close(g);
end

% plot the other things
g=figure; hold on;
bin_centers = bin_edges(1:end-1)+seconds(bin_size/2);
axlabel=[];
for c = 1:length(psth_all)
    if num_trials_all(c) > 15
        plot(bin_centers, psth_all{c})
        axlabel(end+1) = c;
    end
end
xlabel('Time (s)') ; ylabel('Firing Rate (hz)');
legend(label_all(axlabel))

if dosave
    u = goodunit(unit);
    spath = fullfile(savepath, [num2str(u) '--' 'depth' num2str(unit_depth(unit)) '--psths.png']);
    saveas(g, spath);
    close(g);
end

%%

end

%% Plot an example unit ** this is using the old script!!
% align by licks? or align to inferred licks (i.e., if stuck closing doors
% idk)

unit = 154; % filter by depth as well?
dosave = 0;

bychoice = 1;
byreward = 0;

plotOtherMouse = 1;

% plot the other mouse's inferred lick? or reward time?

savepath = 'D:\ANALYSIS\Kevin\Plots\KM33-34\250111\units_by_reward_othermouselick\';
if ~exist(savepath); mkdir(savepath); end

for unit = 1:length(goodunit)

% CHOOSE WHAT TO ALIGN TO

spiketimes = milliseconds(spike_times_byunit{unit});
%trialStarts = milliseconds([info.choice_time{:}]);
trialStarts = milliseconds(choice_time_1); % <- when the choice detectino occurs
trialStarts = milliseconds(choice_time_1_inferred); % <- the first lick
%trialStarts = milliseconds(choice_time_2);

spiketimes_all = []; triallabels_all = [];
for trial = 1:length(trialStarts)
    suse = spiketimes' - trialStarts(trial);
    suse(abs(suse)>seconds(60)) = []; % lets not save every spike every time
    spiketimes_all = [spiketimes_all, suse];
    triallabels_all = [triallabels_all, trial * ones(1,length(suse))];
    %spiketimes_all = [spiketimes_all, spiketimes' - trialStarts(trial)];
    %triallabels_all = [triallabels_all, trial * ones(1,length(spiketimes))];
end

% get other mouse?
if plotOtherMouse
    trialStartOther = milliseconds(choice_time_2_inferred); 
    trialStartOther = trialStartOther - trialStarts;
end


g = figure; hold on;
s = spikeRasterPlot(spiketimes_all,triallabels_all);
s.XLimits = seconds([-5, 10]);

if plotOtherMouse
    hold on;
    plot(1:length(trialStarts), trialStartOther,'ro');
end

title(num2str(unit_depth(unit)))

if dosave
    u = goodunit(unit);
    spath = fullfile(savepath, [num2str(u) '--' 'depth' num2str(unit_depth(unit)) '--alltrials.png']);
    saveas(s, spath);
    close(g);
end

%% NOW PLOT BY GROUPS
% choice and reward and social choice and social reward???
% - maybe condition on one sie?
% otherwise instead of 2x2 im getting a 4x4

% pull in different way so that trials line up for each...
% - and only take trials where BOTH mice made a choice

% load data
LTProb = info.LTProb;
choice = info.choice; choice_time = info.choice_time;
reward = info.reward;
badid = cellfun(@isempty, choice) | cellfun(@isempty, reward) | ...
    cellfun(@length, choice) ~= cellfun(@length, reward);
LTProb(badid) = []; choice(badid) = []; reward(badid) = [];
mouseID = info.mouseID; mouseID(badid) = [];
choice_time(badid) = [];
if info.protocolNum==3
    LT2Prob = info.LT2Prob;
    LT2Prob(badid) = [];
end

% format data
choice1 = zeros(1,length(choice));  % left choice is 1, right is -1!
reward1 = zeros(1,length(choice)); 
choice2 = zeros(1,length(choice)); 
reward2 = zeros(1,length(choice)); 
choice_time1 = zeros(1,length(choice));
choice_time2 = zeros(1,length(choice));
outcomevals = [-1, 1]; % -1 = right, 1 = left
for trial = 1:length(choice)
    for i = 1:length(choice{trial})
        if mouseID{trial}(i) == 0
            % set so assings 1 or -1, and unsassigned is auto 0?
            choice1(trial) = outcomevals( (double(choice{trial}(i)) > 1)+1 );
            reward1(trial) = outcomevals( reward{trial}(i)+1 );
            choice_time1(trial) = choice_time{trial}(i);
        elseif mouseID{trial}(i) == 1
            choice2(trial) = outcomevals( (double(choice{trial}(i)) > 1)+1 );
            reward2(trial) = outcomevals( reward{trial}(i)+1 );
            choice_time2(trial) = choice_time{trial}(i);
        end
    end
end
% filter trials with one choice
badid = choice1==0 | choice2==0;
choice1(badid) = []; choice2(badid) = [];
reward1(badid) = []; reward2(badid) = [];
choice_time1(badid) = []; choice_time2(badid) = [];

% count unique trials to work with
a = [choice1; choice2; reward1; reward2];
[b, ia, ic] = unique(a', 'rows','stable');
h = accumarray(ic, 1);

% maybe just look at choices?
a = [choice1; choice2];
[b, ia, ic] = unique(a', 'rows','stable');
h = accumarray(ic, 1);


% lets do some categories?
% maybe only choice selective...

g = figure;

% Right choice, Right choice
trialID = find(choice1==-1 & choice2==-1);
if byreward; trialID = find(reward1==-1 & reward2==-1); end
subplot(2,2,1);
trial_subset = ismember(triallabels_all, trialID);
RR_spiketimes = spiketimes_all(trial_subset);
RR_spike_trials = triallabels_all(trial_subset);
s = spikeRasterPlot(RR_spiketimes,RR_spike_trials);
s.XLimits = seconds([-5, 10]);
title('Both right')
if byreward; title('Both unrewarded'); end

% self right, other left
trialID = find(choice1==-1 & choice2==1);
if byreward; trialID = find(reward1==-1 & reward2==1); end
subplot(2,2,2);
trial_subset = ismember(triallabels_all, trialID);
RL_spiketimes = spiketimes_all(trial_subset);
RL_spike_trials = triallabels_all(trial_subset);
s = spikeRasterPlot(RL_spiketimes,RL_spike_trials);
s.XLimits = seconds([-5, 10]);
title('Mouse 1 right, mouse 2 left')
if byreward; title('Mouse 1 unrewarded, mouse 2 reward'); end

% both left
trialID = find(choice1==1 & choice2==1);
if byreward; trialID = find(reward1==1 & reward2==-1); end
subplot(2,2,3);
trial_subset = ismember(triallabels_all, trialID);
LL_spiketimes = spiketimes_all(trial_subset);
LL_spike_trials = triallabels_all(trial_subset);
s = spikeRasterPlot(LL_spiketimes,LL_spike_trials);
s.XLimits = seconds([-5, 10]);
title('Both left')
if byreward; title('Mouse 1 reward, mouse 2 unreward'); end

% self left, other right
trialID = find(choice1==1 & choice2==-1);
if byreward; trialID = find(reward1==1 & reward2==1); end
subplot(2,2,4);
trial_subset = ismember(triallabels_all, trialID);
LR_spiketimes = spiketimes_all(trial_subset);
LR_spike_trials = triallabels_all(trial_subset);
s = spikeRasterPlot(LR_spiketimes,LR_spike_trials);
s.XLimits = seconds([-5, 10]);
title('Mouse 1 left, mouse 2 left')
if byreward; title('Both rewarded'); end

if dosave
    spath = fullfile(savepath, [num2str(u) '--' 'depth' num2str(unit_depth(unit)) '-bychoice.png']);
    saveas(g, spath);
    close(g);
end


%% plot firing rates also

nsamps_pergroup = [];


% Right choice, Right choice
trialID = find(choice1==-1 & choice2==-1);
if byreward; trialID = find(reward1==-1 & reward2==-1); end
trial_subset = ismember(triallabels_all, trialID);
% generate psth
bin_size = 0.1;
bin_edges = seconds(-5:bin_size:10);
spike_counts = histcounts(spiketimes_all(trial_subset), bin_edges);
num_trials = length(trialID);
RR_firing_rates = spike_counts / (num_trials * bin_size); 
nsamps_pergroup(1) = num_trials;

% Right choice self, left choice other
trialID = find(choice1==-1 & choice2==1);
if byreward; trialID = find(reward1==-1 & reward2==1); end 
trial_subset = ismember(triallabels_all, trialID);
% generate psth
bin_size = 0.1;
bin_edges = seconds(-5:bin_size:10);
spike_counts = histcounts(spiketimes_all(trial_subset), bin_edges);
num_trials = length(trialID);
RL_firing_rates = spike_counts / (num_trials * bin_size); 
nsamps_pergroup(2) = num_trials;


% Left choice self, right choice other
trialID = find(choice1==1 & choice2==-1);
if byreward; trialID = find(reward1==1 & reward2==-1); end
trial_subset = ismember(triallabels_all, trialID);
% generate psth
bin_size = 0.1;
bin_edges = seconds(-5:bin_size:10);
spike_counts = histcounts(spiketimes_all(trial_subset), bin_edges);
num_trials = length(trialID);
LR_firing_rates = spike_counts / (num_trials * bin_size); 
nsamps_pergroup(3) = num_trials;


% Left choice both
trialID = find(choice1==1 & choice2==1);
if byreward; trialID = find(reward1==1 & reward2==1); end
trial_subset = ismember(triallabels_all, trialID);
% generate psth
bin_size = 0.1;
bin_edges = seconds(-5:bin_size:10);
spike_counts = histcounts(spiketimes_all(trial_subset), bin_edges);
num_trials = length(trialID);
LL_firing_rates = spike_counts / (num_trials * bin_size); 
nsamps_pergroup(4) = num_trials;


g=figure; hold on;
bin_centers = bin_edges(1:end-1)+seconds(bin_size/2);
plot(bin_centers, RR_firing_rates)
plot(bin_centers, RL_firing_rates)
plot(bin_centers, LR_firing_rates)
plot(bin_centers, LL_firing_rates)
xlabel('Time (s)') ; ylabel('Firing Rate (hz)');
legend('Both right','Self right, other left','Self left, other right','Both left')
if byreward
    legend('Neither reward','Mouse 1 unreward, mouse 2 reward','Mouse 1 reward, mouse 2 unreward','Both reward');
end
if dosave
    spath = fullfile(savepath, [num2str(u) '--' 'depth' num2str(unit_depth(unit)) '--firingrates.png']);
    saveas(g, spath);
    close(g);
end

%% Lets split by everything...
g = figure;

count = 1;
for c1 = [-1,1]
    for c2 = [-1,1]
        for r1 = [-1,1]
            for r2 = [-1,1]
                trialID = find(choice1==c1 & choice2==c2 & reward1==r1 & reward2==r2);
                subplot(4,4,count);
                trial_subset = ismember(triallabels_all, trialID);
                RR_spiketimes = spiketimes_all(trial_subset);
                RR_spike_trials = triallabels_all(trial_subset);
                s = spikeRasterPlot(RR_spiketimes,RR_spike_trials);
                s.XLimits = seconds([-5, 10]);
                title(sprintf('choice1: %d, choice2: %d, reward1: %d, reward2: %d', c1,c2,r1,r2));
                count = count+1;
            end
        end
    end
end

%% OK now lets run a significance test
% compare some difference in FR to the distributional difference in FR
% after shuffling 1000 times
% IDK if this is right....

title_label = {'Both right','Self right, other left','Self left, other right','Both left'};
if byreward
    title_label = {'Neither reward','Mouse 1 unreward, mouse 2 reward','Mouse 1 reward, mouse 2 unreward','Both reward'};
end

% compute time range want to analyze
trange = seconds([-4, 4]); % seconds
trangeID = find(bin_centers > trange(1) & bin_centers < trange(2));

% compute average firing rate
spike_counts = histcounts(spiketimes_all, bin_edges);
mu_firing_rates = spike_counts / (length(choice1) * bin_size);
mu_firing_rates = mu_firing_rates(trangeID);

% compute selectivity
group_firing_rates = [RR_firing_rates(trangeID); RL_firing_rates(trangeID); ...
    LR_firing_rates(trangeID); LL_firing_rates(trangeID)];
group_selectivity = max(abs(group_firing_rates' - mu_firing_rates'));

% Now do shuffled!
nedges = [0, cumsum(nsamps_pergroup)];
shuffled_selectivity = [];
for shuff = 1:1000
    % shuffle trials
    trialShuffle = randperm(length(choice1));
    % compute firing rates
    shuffledFR = [];
    for groups = 1:4
        % compute firing rate
        trialID = trialShuffle( (nedges(groups)+1) : (nedges(groups+1)) );
        trial_subset = ismember(triallabels_all, trialID);
        spike_counts = histcounts(spiketimes_all(trial_subset), bin_edges);
        shuffledFR(groups,:) = spike_counts / (length(trialID) * bin_size); 
    end
    % compute selectivity
    shuffledFR = shuffledFR(:, trangeID);
    shuffled_selectivity(shuff, :) = max(abs( shuffledFR' - mu_firing_rates' ));
end

% test plot
g = figure; hold on;
for groups = 1:4
    subplot(2,2,groups); hold on;
    histogram(shuffled_selectivity(:,groups));
    plot(group_selectivity(groups)*[1,1], ylim, 'k--');
    pval = 1 - sum(group_selectivity(groups) > shuffled_selectivity(:,groups)) / length(shuffled_selectivity);
    title(title_label{groups});
    legend(num2str(pval));
end
if dosave
    spath = fullfile(savepath, [num2str(u) '--' 'depth' num2str(unit_depth(unit)) '--selectivity.png']);
    saveas(g, spath);
    close(g);
end


end
