%% Regression for choice, reward, value
% including data from multipel recordings and multiple sessions

% 1) choice and reward and interaction
% - balance classes!
% 1.5) social as well
% - do I want to compare the social and non social fit strengths?

% 2) Value
% - idk if need to balance classes?

% 3) Shuffled control
% - hard to do 'in place', where keep task statistics the same
% - since have a limited number of swaps per session
% - maybe just swap the choice labels randomly, rather than in place

%% Paths

basepath = 'Z:\HarveyLab\Tier1\Kevin\Videos\';

% miceNonsocial = {'KM33','KM35','KM37'};
% dates_nonsocial = {{'250207','250213','250221'}; ...
%     {'250226','250307','250314'}; ...
%     {'250227','250306','250313'}};
% 
% miceSocial = {'KM33-34','KM35-36','KM37-38'};
% dates_social = {{'250208','250214','250222'}; ...
%     {'250227','250309','250315'}; ...
%     {'250228','250307','250314'}};

% social plus
miceUse= {'KM33-34','KM35-36','KM37-38'};
dateUse = {{'250208','250214','250222'}; ...
    {'250227','250309','250315'}; ...
    {'250228','250307','250314'}};

% non social version of above?
miceUse= {'KM33','KM35','KM37'};
dateUse = {{'250207','250213','250221'}; ...
    {'250226','250307','250314'}; ...
    {'250227','250306','250313'}};

% social minus (to do!)
miceUse= {'KM33-34','KM35-36','KM37-38'};
dateUse = {{'250307','250313','250320'};...% chosen 1 week apart
    {'250328','250403','250409'};...% first one might be early?
    {'250325','250331','250405'}};

% 
% social neutral
miceUse = {'KM33-34', 'KM35-36','KM37-38'};
dateUse = {{'250418'},...
    {'250430'},{'250417'}}; 

dateUse = {{'250408','250418'};...
    {'250416','250430'};...
    {'250408','250417'}};

% New regions
miceUse = {'KM41-42'};
dateUse = {{'250420','250429','250509'}};

% VCA1 only one mouse for now. maybe average over two sesion days?
miceUse = {'KM41-42','KM41-42'};
dateUse = {{'250730'},{'250805'}};

% Lets try these other regions again, but for more mice
%-  although I think I missed a bunch
miceUse = {'KM41-42','KM43-44','KM45-46'};
dateUse = {{'250429'};...
    {'250731'};...
    {'250802'}};

% Non social version of VCA1 for sanity?
miceUse = {'KM42','KM42'};
dateUse = {{'250729'},{'250804'}};

miceUse = {'KM47-48'};
dateUse = {{'250809'}};

% KM45 only in BLA?
miceUse = {'KM45-46'};
dateUse = {{'250820','250904'}};

%% ** NEW GROUP BY REGIONS:
% - this does for all mice, but can filter by depht later?
% KM41, in PPC-CA2-BLA
% KM42, in VCA1
% KM43, in PPC-CA2-BLA
% KM44, nothing
% KM45, in BLA
% KM46, in ?? (probably PPC-CA2-BLA)
% KM47, in VCA1
% KM48, in VCA1 and LEC

% 0) mPFC!
miceUse = {'KM33-34','KM35-36','KM37-38'};
dateUse = {{'250208','250214','250222'}; ...
    {'250227','250309','250315'}; ...
    {'250228','250307','250314'}};
mouseIDuse = [1, 1, 1];
depth_range  = {[640, 2340], [570, 2410], [660, 2430]}; % range on NP in mm, relative to top of brain
top_of_brain = [383, 383, 370]; % <- which probe site is at the top of the brain
active_probe_range_Allen = [489, 436, 490];

%** Everything below is a guess for Chris
% 1) BLA
miceUse = {'KM41-42','KM43-44','KM45-46','KM45-46'};
dateUse = {{'250512','250505'}; ...
    {'250819','250806','250731'}; ...
    {'250909','250820'}; ...
    {'250906','250902','250813'}}; % <- hmm
mouseIDuse = [1, 1, 1, 2]; % which mouse in pair was recored from
depth_range = {[5310, 6050], [4600, 6000], [4600, 6000], [4600, 6000],}; % <- this is borders relative to top!!! <- s
% somethings up with the top of brain and the allen range...
top_of_brain = [562, 383+114, 383+113, 383+223];
active_probe_range_Allen = [624, 383+114, 383+113, 383+223]; % no scaling

% - need to coordinate relative to 

% 2) VCA1
miceUse = {'KM41-42', 'KM47-48','KM47-48'};
dateUse = {{'250728','250805'};...
    {'250821','250911','250903'};...
    {'250909','250915'}};
mouseIDuse = [2, 1, 2];
depth_range = {[4380, 4740], [3200, 3830],[3200, 3830]}; % 
top_of_brain = [383+96, 383, 383+19];
active_probe_range_Allen = [500, 383, 383+19]; % no scaling

% 3) CA2
miceUse = {'KM41-42','KM43-44', 'KM45-46'};
mouseIDuse = [1, 1, 2]
depth_range = {[1200, 2500], [1200, 2500], [1200, 2500]}; % <- this is borders relative to top!!!
active_probe_range_Allen = top_of_brain; % no scaling

% 4) PPC?
miceUse = {'KM41-42','KM43-44', 'KM45-46'};
mouseIDuse = [1, 1, 2]
depth_range = {[0, 1000], [0, 1000], [0, 1000]}; % <- this is borders relative to top!!!
active_probe_range_Allen = top_of_brain; % no scaling

% IF want to include everything
% just make depth_range -inf to inf

%% Paralellize
% miceUse= {'KM33-34','KM33-34','KM33-34','KM35-36','KM35-36','KM35-36'...
%     ,'KM37-38','KM37-38','KM37-38'};
% dateUse = {'250208','250214','250222', ...
%     '250227','250309','250315', ...
%     '250228','250307','250314'};

% straightforward test
% session_paths = 'Z:\HarveyLab\Tier1\Kevin\Videos\KM33-34\250222_g0';
% this was run on other script (test_regression_analysis.m)
% so lets try on this one?
basepath = 'Z:\HarveyLab\Tier1\Kevin\Videos\';
addpath(genpath('Z:\HarveyLab\Tier1\Kevin\Analysis\20250718_backup_Cindys_PC\Utilities'))
addpath(genpath('Z:\HarveyLab\Tier1\Kevin\Analysis\20250718_backup_Cindys_PC\RL_modeling'))
miceUse = {'KM33-34'};
dateUse = {{'250222'}};

miceUse = {'KM41-42'};
dateUse = {{'250805'}};
%% Variables to save


% Saved variables
% Choice and reward
p_mdl = []; % time x term x neuron
p_mdl_social_unbalanced = [];
p_mdl_social_balanced_choice = [];
p_mdl_social_balanced_reward = [];
b_mdl = []; % time x term x neuron
b_mdl_social_unbalanced = [];
b_mdl_social_balanced_choice = [];
b_mdl_social_balanced_reward = [];

p_mdl_social_shuffled = [];
b_mdl_social_shuffled = [];
p_mdl_social_shuffled_control = []; 
b_mdl_social_shuffled_control = [];

% samples in balanced
minN_balance_save = [];

% history
p_mdl_hist = [];
b_mdl_hist = [];

p_mdl_hist_2 = []; % two back for a example fig
b_mdl_hist_2 = []; 

% Value
p_mdl_Q = [];
p_mdl_Q_social = [];
b_mdl_Q = [];
b_mdl_Q_social = [];

p_mdl_Q_social_shuffled = [];
p_mdl_Q_social_shuffled_control = [];
b_mdl_Q_social_shuffled = [];
b_mdl_Q_social_shuffled_control = [];

% joint
p_mdl_joint = [];
b_mdl_joint = [];
p_mdl_social_joint = [];
b_mdl_social_joint = [];

% interaction??
p_mdl_interaction = [];
b_mdl_interaction = [];
p_mdl_interaction_social = [];
b_mdl_interaction_social = [];



mouseIDsave = [];
dateIDsave = [];
unitIDsave = [];
mouseNamesave = {};
dateNamesave = {};

unitdepth = [];

% want to save all the activity?

%% settings

bin_size = 0.2; % 200 ms bins
bin_edges = (-3:bin_size:5); % time around lick
nShuff = 10; % number of times to rebalance classes randomly

do_glm = 0; % flag to do linear model of glm

do_zscore_activity = 0; % flag to normalize activity (over entire bin_edge range)
norm_window = [-10, bin_edges(1)]; % ?

model_use = 'non-social';

use_other_mouses_choicetime_as_0 = 0;

goodid_criterion = 1; % 1 is the strict combo, 0 is the less strict allen one?


which_mouse_has_implant = 2; % or 2 for the pair of the mouse!
% *** 


%% Loop

% counters
ucount = 1;

% loop
for m = 1:length(miceUse)
    for d = 1:length(dateUse{m})

session_path = fullfile(basepath, miceUse{m}, [dateUse{m}{d} '_g0']);
%session_path = fullfile(basepath, miceUse{m}, [dateUse{m} '_g0']);

disp(session_path)
%% Load spike info

% extract mouse name
mouse_name = strsplit(session_path, '\');
mouse_name = mouse_name(contains(mouse_name, 'KM'));

% load spike times
spike_times = readNPY(fullfile(session_path, 'Spike_Sorting','phy', 'spike_times.npy'));
spike_clusters = readNPY(fullfile(session_path, 'Spike_Sorting', 'phy', 'spike_clusters.npy'));
T_1 = parquetread(fullfile(session_path, 'Spike_Sorting', 'quality_metrics.parquet'));
T_1 = renamevars(T_1, 'firing_rate', 'fr');
channel_positions = readNPY(fullfile(session_path,'Spike_Sorting', 'phy', 'channel_positions.npy'));
primary_channels = readNPY(fullfile(session_path,'Spike_Sorting', 'primary_channels.npy'));
depth_1 = channel_positions(primary_channels+1 , 2);

% load TTL and sync
load(fullfile(session_path, 'TTLs','TTLs.mat'));
spike_times_NIDQ = NPtoNIDQ(1)*double(spike_times) + NPtoNIDQ(2);
spike_times_teensy = NIDQtoTeensy(1)*spike_times_NIDQ + NIDQtoTeensy(2);

% reformat data per unit (watch the index by zero!!)
%goodid = T_1.snr>4.0 & T_1.fr>0.01 & T_1.isi_violations_ratio<1; % <- isabels original
%goodid = T_1.isi_violations_ratio<1 & T_1.amplitude_cutoff<0.1 & T_1.presence_ratio>0.9; % Allen, and Isabel also uses
if ~goodid_criterion
goodid = T_1.snr>1.0 & T_1.fr>0.01 & T_1.isi_violations_ratio<1 & T_1.presence_ratio>0.5; % allen with some actual fr thresholds
else
goodid = (T_1.snr > 2.0 & T_1.fr > 0.05 & T_1.nn_hit_rate > 0.5 & ...
           T_1.isi_violations_ratio<1 & T_1.amplitude_cutoff<0.1 & T_1.presence_ratio>0.9); % Allen + some extra from isabel
end
disp(['num units passing criterion: ' num2str(sum(goodid))])

goodid = find(goodid);

spike_times_byunit = cell(max(spike_clusters)+1,1);
for i = 0:max(spike_clusters)
    spike_times_byunit{i+1} = spike_times_teensy(spike_clusters==i);
end


% load behavioral data
beh_path = [session_path(1:end-3) '*.dat'];
if contains(mouse_name, '-')
    info = createBEHstruct_social(mouse_name{1}, dir( beh_path));
else
    info = createBEHstruct_nonsocial(mouse_name{1}, dir(beh_path));
end

% %% TODO: filter by depth to separate brain regions
% % Using *only* SharpTrack, get the coordinates along the probe that we care
% % about
% 
% % Convert probe points from tip to brain surface
% depth_fromtop = depth_1(goodid) - top_of_brain(m)*10;
% % Scale from Allen CCF to emperical probe units
% % - probe_range is how deep the tip goes based on final point. 
% depth_in_AllenCCF = depth_fromtop * (active_probe_range_Allen(m) / top_of_brain(m)) * -1; 
% % Keep those in region of interest along probe
% ROI = find(depth_in_AllenCCF < depth_range{m}(2) & depth_in_AllenCCF > depth_range{m}(1));
% goodid = goodid(ROI);
% 
% % Set up which mouse of pair to use for aligment
% which_mouse_has_implant = mouseIDuse(m); % new
% 
% 
% disp(['num units in correct brain area: ' num2str(length(goodid))])

%% Alternative: just use a heuristic based on what I see along probe
% e.g., its p easy to spot hippocampus


%% 1) First encoding model
% with balanced classes

% pull behavior information
if contains(mouse_name, '-')
    if which_mouse_has_implant==1
        [choice, reward, choice_time, choice2, reward2, choice_time2, LTProb] = extract_session_params(info);
    elseif which_mouse_has_implant==2
        [choice2, reward2,choice_time2, choice,reward,choice_time, LTProb] = extract_session_params(info);
    else
        disp('paramter wrong!')
    end
else
[choice, reward, choice_time] = extract_session_params(info);
end
% pull value information!
[Q,p, acc] = extract_value_information(info, model_use);
%[Q,p, acc] = extract_value_information(info, model_use, 1, 1, 1);

if contains(mouse_name, '-')
    if which_mouse_has_implant == 1
        [Q2, ~, acc] = extract_value_information(info, model_use, 2); 
        %[Q2, ~, acc] = extract_value_information(info, model_use, 2, 1, 1); 
    elseif which_mouse_has_implant == 2
        [Q, ~, acc] = extract_value_information(info, model_use, 2); 
        [Q2,p, acc] = extract_value_information(info, model_use);

    else
        disp('parameter wrong!')
    end
end


% filter for crashed recordings (but this should be fine for any)
crash_thresh = max(spike_times_teensy);
if any(choice_time > crash_thresh)
    % trials were recorded after the last found spike!
    disp('there was a crash in this recording!')
    crash_trial = min(find(choice_time > crash_thresh));
    choice(:,crash_trial:end) = [];
    reward(:,crash_trial:end) = [];
    choice_time(:,crash_trial:end) = [];
    Q(:,crash_trial:end) = [];
    LTProb(:,crash_trial:end) = [];
    if contains(mouse_name, '-')
        choice2(:,crash_trial:end) = [];
        reward2(:,crash_trial:end) = [];
        choice_time2(:,crash_trial:end) = [];
        Q2(:,crash_trial:end) = [];
    end
end


% loop over units
for u = 1:length(goodid)

    unit = goodid(u); 
    disp(u/length(goodid))

    % pull spike times
    spiketimes = milliseconds(spike_times_byunit{unit});
    trialStarts = milliseconds(choice_time);

    if use_other_mouses_choicetime_as_0
        trialStarts = milliseconds(choice_time2);
    end
    
    spiketimes_all = {}; % Mx1 cell array of spike times
    for trial = 1:length(trialStarts)
        suse = spiketimes' - trialStarts(trial);
        suse(abs(suse)>seconds(60)) = []; % lets not save every spike every time
        spiketimes_all{trial} = seconds(  suse  );
    end

    % loop over each time bin in range
    for j = 1:(length(bin_edges)-1)
        % get activity(t)
        spike_counts = cellfun(@(v) sum(v >= bin_edges(j) & v < bin_edges(j+1)), spiketimes_all);
        activity = spike_counts ./ bin_size;

        % z-score activity using entire trial window (maybe use even larger
        % range?)
        if do_zscore_activity
            %spike_norm = cellfun(@(v) sum(v >= bin_edges(1) & v < bin_edges(end)), spiketimes_all);
            %activity_baseline = spike_norm ./ (bin_edges(end) - bin_edges(1));
            spike_norm = cellfun(@(v) sum(v >= norm_window(1) & v < norm_window(end)), spiketimes_all);
            activity_baseline = spike_norm ./ (norm_window(end) - norm_window(1));
            mu_activity = mean(activity_baseline);
            std_activity = std(activity_baseline);
            activity = (activity - mu_activity) / std_activity;
        end
        
        % 1) chocie encoding
        % fit multiple linear regression 
        X = [choice', reward', ]; % include history model?
        if do_glm
            mdl = fitglm(X, activity, 'Distribution','poisson');
        else
            mdl = fitlm(X, activity);
        end
        % save
        p_mdl(j, :, ucount) = mdl.Coefficients.pValue(2:end);
        b_mdl(j, :, ucount) = mdl.Coefficients.Estimate(2:end);

        % 1.5) choice history?
        Xhist = [X(2:end,:), X(1:end-1, :)];
        if do_glm
            mdl_hist = fitglm(Xhist, activity(2:end),'Distribution','poisson');
        else
            mdl_hist = fitlm(Xhist, activity(2:end));
        end
        p_mdl_hist(j,:,ucount) = mdl_hist.Coefficients.pValue(2:end);
        b_mdl_hist(j,:,ucount) = mdl_hist.Coefficients.Estimate(2:end);

        % 1.75) choice history 2x
        Xhist2 = [X(3:end,:), X(2:end-1, :), X(1:end-2,:)];
        if do_glm
            mdl_hist2 = fitglm(Xhist2, activity(3:end), 'Distribution','poisson');
        else
            mdl_hist2 = fitlm(Xhist2, activity(3:end));
        end
        p_mdl_hist_2(j,:,ucount) = mdl_hist2.Coefficients.pValue(2:end);
        b_mdl_hist_2(j,:,ucount) = mdl_hist2.Coefficients.Estimate(2:end);

        % regression for the social mice
        if contains(mouse_name, '-')
            Xsocial = [choice', reward', choice2', reward2'];
            if do_glm
                mdl_social = fitglm(Xsocial, activity,'Distribution','Poisson');
            else
                mdl_social = fitlm(Xsocial, activity);
            end

            p_mdl_social_unbalanced(j,:,ucount) = mdl_social.Coefficients.pValue(2:end);
            b_mdl_social_unbalanced(j,:,ucount) = mdl_social.Coefficients.Estimate(2:end);
            
            % Shuffle 1: balanced classes
            % - pros - removes colinearity
            % - cons - not a lot of samples left, which affects encoding 
            for n = 1:nShuff

                % subselect to balance choice classes
                [uniqueChoices, ~, choiceIDX] = unique([choice; choice2]', 'rows');
                minN = min(accumarray(choiceIDX,1));
                subselectIdx = [];
                for i = 1:size(uniqueChoices,1)
                    trialIndices = find(choiceIDX==i);
                    subselectIdx = [subselectIdx; randsample(trialIndices, minN)];
                end
                
                if do_glm
                    mdl_social = fitglm(Xsocial(subselectIdx, :), activity(subselectIdx),'Distribution','poisson');
                else
                    mdl_social = fitlm(Xsocial(subselectIdx, :), activity(subselectIdx));
                end
                p_mdl_social_balanced_choice(j,:,ucount,n) = mdl_social.Coefficients.pValue(2:end);
                b_mdl_social_balanced_choice(j,:,ucount,n) = mdl_social.Coefficients.Estimate(2:end);

                % subselect to alance reward
                [uniqueRewards, ~, rewardIDX] = unique([reward; reward2]', 'rows');
                minN = min(accumarray(rewardIDX, 1));
                subselectIdx = [];
                for i = 1:length(uniqueRewards)
                    trialIndices = find(rewardIDX==i); % is this right?
                    subselectIdx = [subselectIdx; randsample(trialIndices, minN)];
                end
                if do_glm
                    mdl_social = fitglm(Xsocial(subselectIdx, :), activity(subselectIdx),'Distribution','poisson');
                else
                    mdl_social = fitlm(Xsocial(subselectIdx, :), activity(subselectIdx));
                end
                p_mdl_social_balanced_reward(j,:,ucount,n) = mdl_social.Coefficients.pValue(2:end);
                b_mdl_social_balanced_reward(j,:,ucount,n) = mdl_social.Coefficients.Estimate(2:end);

            end % shuffle

            % Shuffle 2: keep the statistics intact
            % - dont need to do even odd splits, just need to flip choice
            % around
    
            % Generate groupings
            minN = min(diff(find(abs(diff(LTProb))>0))) - 1; % yeah yeah this is terrible i agree
            minN = min(minN, 30); % um. 
            swapid = find(abs(diff(LTProb))>20);
            swapparity = 1:length(swapid);
            % toss swaps out of range
            badid = find(diff(swapid) < minN); % swaps that arent long enough
            badid = [badid, find(swapid - minN <= 0), find(swapid + minN > length(choice))]; % swaps at start and end without enough trials
            swapid(badid) = [];
            swapparity(badid) = [];
            for n = 1:nShuff
                % swap and generate data
                shuff = randperm(length(swapid));
                shuffid = swapid(shuff);
                shuffparity = swapparity(shuff);
                Xswap = []; 
                activity_swap = [];
                for sid = 1:length(swapid)
                    qq = swapid(sid) + (-minN : minN);
                    rr = shuffid(sid) + (-minN : minN);
                    % if even, flip choice parity. if odd, do nothing
                    parity = 2*mod(shuffparity(sid),2) - 1;
                    Xswap = [Xswap; choice(qq)', reward(qq)', choice2(rr)'*parity, reward2(rr)'];
                    activity_swap = [activity_swap, activity(qq)];
                end
                if do_glm
                mdl_shuff = fitglm(Xswap, activity_swap,'Distribution','poisson');
                else
                mdl_shuff = fitlm(Xswap, activity_swap);
                end
                p_mdl_social_shuffled(j,:,ucount, n) = mdl_shuff.Coefficients.pValue(2:end);
                b_mdl_social_shuffled(j,:,ucount, n) = mdl_shuff.Coefficients.Estimate(2:end);
            end % end shuffle
            % run one model grouped the same way as a control
            Xswap = []; 
            activity_swap = [];
            for sid = 1:length(swapid)
                qq = swapid(sid) + (-minN : minN);
                Xswap = [Xswap; choice(qq)', reward(qq)', choice2(qq)', reward2(qq)'];
                activity_swap = [activity_swap, activity(qq)];
            end
            if do_glm
                mdl_shuff = fitglm(Xswap, activity_swap,'Distribution','poisson');
            else
                mdl_shuff = fitlm(Xswap, activity_swap);
            end
            p_mdl_social_shuffled_control(j,:,ucount) = mdl_shuff.Coefficients.pValue(2:end);
            b_mdl_social_shuffled_control(j,:,ucount) = mdl_shuff.Coefficients.Estimate(2:end);

        end % social if statement

        % 2) value encoding
        deltaQ = Q(1,:) - Q(2,:); % contra - ipsilateral (right - left)
        X = [deltaQ'];
        if do_glm
            mdl = fitglm(X, activity, 'Distribution','poisson');
        else
            mdl = fitlm(X, activity);
        end
        p_mdl_Q(j, :, ucount) = mdl.Coefficients.pValue(2);
        b_mdl_Q(j, :, ucount) = mdl.Coefficients.Estimate(2);

%         % joint value encoding!!!
%         X = [choice', reward', deltaQ'];
%         if do_glm
%             mdl = fitglm(X, activity, 'Distribution','poisson');
%         else
%             mdl = fitlm(X, activity);
%         end
%         p_mdl_joint(j, :, ucount) = mdl.Coefficients.pValue(2);
%         b_mdl_joint(j, :, ucount) = mdl.Coefficients.Estimate(2);

        % social
        if contains(mouse_name, '-')
            deltaQ2 = Q2(1,:) - Q2(2,:);
            X = [deltaQ', deltaQ2'];
            if do_glm
                mdl = fitglm(X, activity, 'Distribution','poisson');
            else
                mdl = fitlm(X, activity);
            end
            p_mdl_Q_social(j,:,ucount) = mdl.Coefficients.pValue(2:end);
            b_mdl_Q_social(j,:,ucount) = mdl.Coefficients.Estimate(2:end);

            % shuffle here as well?
            % - do the statistics shuffle... idk
            % Generate groupings
            minN = min(diff(find(abs(diff(LTProb))>0))) - 1; % yeah yeah this is terrible i agree
            minN = min(minN, 30); % um. 
            swapid = find(abs(diff(LTProb))>20);
            swapparity = 1:length(swapid);
            % toss swaps out of range
            badid = find(diff(swapid) < minN); % swaps that arent long enough
            badid = [badid, find(swapid - minN <= 0), find(swapid + minN > length(choice))]; % swaps at start and end without enough trials
            swapid(badid) = [];
            swapparity(badid) = [];
            for n = 1:nShuff
                % swap and generate data
                shuff = randperm(length(swapid));
                shuffid = swapid(shuff);
                shuffparity = swapparity(shuff);
                Xswap = []; 
                activity_swap = [];
                for sid = 1:length(swapid)
                    qq = swapid(sid) + (-minN : minN);
                    rr = shuffid(sid) + (-minN : minN);
                    % if even, flip choice parity. if odd, do nothing
                    parity = 2*mod(shuffparity(sid),2) - 1;
                    Xswap = [Xswap; deltaQ(qq)', deltaQ2(rr)'*parity]; %<- this isnt right...
                    activity_swap = [activity_swap, activity(qq)];
                end
                if do_glm
                mdl_shuff = fitglm(Xswap, activity_swap,'Distribution','poisson');
                else
                mdl_shuff = fitlm(Xswap, activity_swap);
                end
                p_mdl_Q_social_shuffled(j,:,ucount, n) = mdl_shuff.Coefficients.pValue(2:end);
                b_mdl_Q_social_shuffled(j,:,ucount, n) = mdl_shuff.Coefficients.Estimate(2:end);
            end % end shuffle
            % run one model grouped the same way as a control
            Xswap = []; 
            activity_swap = [];
            for sid = 1:length(swapid)
                qq = swapid(sid) + (-minN : minN);
                Xswap = [Xswap; deltaQ(qq)', deltaQ2(qq)'];
                activity_swap = [activity_swap, activity(qq)];
            end
            if do_glm
                mdl_shuff = fitglm(Xswap, activity_swap,'Distribution','poisson');
            else
                mdl_shuff = fitlm(Xswap, activity_swap);
            end
            p_mdl_Q_social_shuffled_control(j,:,ucount) = mdl_shuff.Coefficients.pValue(2:end);
            b_mdl_Q_social_shuffled_control(j,:,ucount) = mdl_shuff.Coefficients.Estimate(2:end);
            

        end

        % 3) joint choice, reward, value encoding
        % - just do the damn unshuffled, unbalanced version
        % - I think this will be useful to have
        % fit multiple linear regression 
        deltaQ = Q(1,:) - Q(2,:); 
        X = [choice', reward', deltaQ' ]; % include history model?
        if do_glm
            mdl = fitglm(X, activity, 'Distribution','poisson');
        else
            mdl = fitlm(X, activity);
        end
        % save
        p_mdl_joint(j, :, ucount) = mdl.Coefficients.pValue(2:end);
        b_mdl_joint(j, :, ucount) = mdl.Coefficients.Estimate(2:end);
        if contains(mouse_name, '-')
            deltaQ2 = Q2(1,:) - Q2(2,:);
            Xsocial = [choice', reward', deltaQ', choice2', reward2', deltaQ2'];
            if do_glm
                mdl_social = fitglm(Xsocial, activity,'Distribution','Poisson');
            else
                mdl_social = fitlm(Xsocial, activity);
            end

            p_mdl_social_joint(j,:,ucount) = mdl_social.Coefficients.pValue(2:end);
            b_mdl_social_joint(j,:,ucount) = mdl_social.Coefficients.Estimate(2:end);
        end

        % 4) interactions encoding models
        % - this includes interaction terms of decision variables (not reward!)
        X = [choice', reward', (choice.*reward)'];
        if do_glm
            mdl_interaction = fitglm(X, activity, 'Distribution','poisson');
        else
            mdl_interaction = fitlm(X, activity);
        end
        p_mdl_interaction(j,:,ucount) = mdl_interaction.Coefficients.pValue(2:end);
        b_mdl_interaction(j,:,ucount) = mdl_interaction.Coefficients.Estimate(2:end);
        if contains(mouse_name, '-')
            X = [choice', reward', (choice.*reward)', choice2', reward2', (choice2.*reward2)',...
                (choice.*choice2)', (reward.*reward2)', (choice.*choice2.*reward.*reward2)' ]; 
            % interaction_labels_social = {'Choice','Reward','Interaction','Social choice','Social reward',...
            % 'Social interaction','Joint choice','Joint reward','Joint everything!'};
            if do_glm
                mdl_interaction_social = fitglm(X, activity,'Distribution','poisson');
            else
                mdl_interaction_social = fitlm(X, activity);
            end
            p_mdl_interaction_social(j,:,ucount) = mdl_interaction_social.Coefficients.pValue(2:end);
            b_mdl_interaction_social(j,:,ucount) = mdl_interaction_social.Coefficients.Estimate(2:end);
        end

        % 5) interactison in values?
        % - coudl these be interpreted as something interesting? i.e.,
        % reward inequality or something
        % Qchoice, Qchoice2, Qchoicedifference;
        % deltaQ1 + deltaQ2; deltaQ1 - deltaQ2;  < - these should show up
        
        
    end % time bin edges

    mouseIDsave(ucount) = m;
    dateIDsave(ucount) = d;
    unitIDsave(ucount) = unit;
    unitdepth(ucount) = depth_1(unit);
    minN_balance_save(ucount) = minN;
    mouseNamesave{ucount} = miceUse{m};
    dateNamesave{ucount} = dateUse{m}{d};
    %dateNamesave{ucount} = dateUse{m};%{d};
    ucount = ucount+1;
end % unit



        


    end % dates
end % mice



%% Save these variables since everythign takes so long to run?
% 
% savepath = 'D:\MEETINGS\Kevin\250527\encoding_data_nonsocial.mat';
 %savepath = 'D:\MEETINGS\Kevin\250527\encoding_data_bla_hc_ppc.mat';
%savepath = 'D:\MEETINGS\Kevin\250527\encoding_data_social_minus_shufflefixed.mat';
savepath = 'D:\MEETINGS\Kevin\250527\encoding_data_social_plus_shufflefixed.mat';
%savepath = 'D:\MEETINGS\Kevin\250527\encoding_data_social_neutral.mat';
%savepath = 'D:\MEETINGS\Kevin\250527\encoding_data_nonsocial_withjoint.mat';
%savepath = 'D:\MEETINGS\Kevin\250527\encoding_data_social_plus_zscored.mat';
%savepath = 'D:\MEETINGS\Kevin\250527\encoding_data_bla_hc_ppc_socialplus.mat';
savepath = 'D:\MEETINGS\Kevin\250527\encoding_data_social_neutral_moredates.mat';

%savepath = 'D:\MEETINGS\Kevin\250527\encoding_data_social_aligned_to_other_mouse.mat';

savepath = 'D:\MEETINGS\Kevin\251004\VCA1_first_pass.mat';
savepath = 'D:\MEETINGS\Kevin\251004\bla_hc_ppc_moremice_maybe.mat';
savepath = 'D:\MEETINGS\Kevin\251004\VCA1_second_pass.mat';
savepath = 'D:\MEETINGS\Kevin\251004\VCA1_km47_first_pass.mat';


savepath = 'D:\MEETINGS\Kevin\251004\BLA_KM45_first_pass.mat';


savepath = 'D:\MEETINGS\Kevin\251004\mPFC_socialplus_with_histology.mat';
savepath = 'D:\MEETINGS\Kevin\251004\BLA_socialplus_makinguphistology_TEMP.mat';
savepath = 'D:\MEETINGS\Kevin\251004\VCA1_socialplus_makinguphistology_TEMP.mat';


savepath = 'D:\MEETINGS\Kevin\251004\PPC-HC-BLA_socialplus_everything_TEMP.mat';


savepath= 'Z:\HarveyLab\Tier1\Bing_Shiuan\Codes\250925_41-42_250806.mat'
if ~exist(savepath)
save(savepath,...
    'mouseIDsave','dateIDsave','unitIDsave','unitdepth',...
    'miceUse','dateUse','basepath',...
    'bin_size','bin_edges','nShuff','do_glm',...
    'p_mdl','p_mdl_social_unbalanced','p_mdl_social_balanced_choice','p_mdl_social_balanced_reward',...
    'b_mdl','b_mdl_social_unbalanced','b_mdl_social_balanced_choice','b_mdl_social_balanced_reward',...
    'p_mdl_hist','b_mdl_hist','p_mdl_Q','b_mdl_Q','p_mdl_Q_social','b_mdl_Q_social',...
    'p_mdl_social_shuffled','p_mdl_social_shuffled_control','b_mdl_social_shuffled','b_mdl_social_shuffled_control',...
    'p_mdl_Q_social_shuffled','p_mdl_Q_social_shuffled_control','b_mdl_Q_social_shuffled','b_mdl_Q_social_shuffled_control',...
    'p_mdl_social_joint','b_mdl_social_joint','p_mdl_joint','b_mdl_joint',...
    'p_mdl_interaction_social','b_mdl_interaction_social','p_mdl_interaction','b_mdl_interaction');
else
    disp('may be overwriting some data!')
end


%%
% save('D:\MEETINGS\Kevin\250527\encoding_data_social.mat',...
%    'p_mdl_interaction_social','b_mdl_interaction_social','p_mdl_interaction','b_mdl_interaction',...
%     '-append')

%% Load other data
load('D:\MEETINGS\Kevin\250527\encoding_data_social_plus_shufflefixed.mat')
load('D:\MEETINGS\Kevin\250527\encoding_data_nonsocial.mat')


%%
%load('D:\MEETINGS\Kevin\250527\encoding_data.mat')
%loadpath = 'D:\MEETINGS\Kevin\250527\encoding_data_bla_hc_ppc_socialplus.mat';
% loadpath = 'D:\MEETINGS\Kevin\251004\bla_hc_ppc_moremice_maybe.mat';
% loadpath = 'D:\MEETINGS\Kevin\251004\PPC-HC-BLA_socialplus_everything_TEMP.mat';

load('Z:\HarveyLab\Tier1\Bing_Shiuan\Codes\250925.mat');

% crop for a brain region based on depth?
unitdepth_fromtop = unitdepth - max(unitdepth);
hc_range = [-1.2, -2.5]*1000;
bla_range = [-4.6, -6]*1000;
ppc_range = [0, -1.0]*1000;


region_use = find(unitdepth_fromtop < hc_range(1) & unitdepth_fromtop > hc_range(2)); 
%region_use = find(unitdepth_fromtop < bla_range(1) & unitdepth_fromtop > bla_range(2)); 
%region_use = find(unitdepth_fromtop < ppc_range(1) & unitdepth_fromtop > ppc_range(2)); 

% keep only neurons want to use for all saved data
% then resave and just load that?
varNames = { 'mouseIDsave','dateIDsave','unitIDsave','unitdepth',...
    'p_mdl','p_mdl_social_unbalanced','p_mdl_social_balanced_choice','p_mdl_social_balanced_reward',...
    'b_mdl','b_mdl_social_unbalanced','b_mdl_social_balanced_choice','b_mdl_social_balanced_reward',...
    'p_mdl_hist','b_mdl_hist','p_mdl_Q','b_mdl_Q','p_mdl_Q_social','b_mdl_Q_social',...
    'p_mdl_social_shuffled','p_mdl_social_shuffled_control','b_mdl_social_shuffled','b_mdl_social_shuffled_control',...
    'p_mdl_Q_social_shuffled','p_mdl_Q_social_shuffled_control','b_mdl_Q_social_shuffled','b_mdl_Q_social_shuffled_control',...
    'p_mdl_social_joint','b_mdl_social_joint','p_mdl_joint','b_mdl_joint',...
    'p_mdl_interaction_social','b_mdl_interaction_social','p_mdl_interaction','b_mdl_interaction'};

nNeurons = length(unitdepth);
for i = 1:length(varNames)
    if exist(varNames{i}, 'var') 
        data = evalin('base', varNames{i}); % Grab the variable
        dataSize = size(data);
        dim_to_crop = find(dataSize == nNeurons);
        if ~isempty(dim_to_crop)
            % Only crop along the first matching dimension
            dim = dim_to_crop(1);
            
            % Build dynamic indexing
            idx = repmat({':'}, 1, ndims(data));
            idx{dim} = region_use;

            % Apply cropping
            data = data(idx{:});
            assignin('base', varNames{i}, data); % Put cropped data back into workspace
        else
            fprintf('Skipping "%s": no dimension matches length of region_use.\n', varNames{i});
        end    
    end
end

%% Need some clean way of doing this for each mouse with histology?


%% Plots
% 1) Non-social plots! of fraction selective neurons

% a) fraction significant neurons over time
labels = {'Choice','Reward'};
frac_sig = zeros(length(bin_edges)-1, 2);
for term = 1:2
    num_significant = squeeze(p_mdl(:, term, :) < 0.05);
    frac_sig(:,term) = sum(num_significant, 2) / length(unitIDsave);
end

bin_centers = (bin_edges(1:end-1) + bin_edges(2:end))/2;
% figure;
% plot(bin_centers, frac_sig); 
% legend(labels);
% xlabel('Time to choice (s)'); ylabel('Fraction selective')

% do by mouse? (and session?)
frac_sig_mouse = zeros(length(bin_edges)-1, 2, length(miceUse));
for term = 1:2
    for m = 1:length(miceUse)
        num_significant = squeeze(p_mdl(:, term, mouseIDsave==m) < 0.05);
        frac_sig_mouse(:, term, m) = sum(num_significant, 2) / sum(mouseIDsave==m);
    end
end

% plot
figure; hold on; coloruse = 'br'; axlb = [];
for term = 1:2
    for m = 1:length(miceUse)
        plot(bin_centers, frac_sig_mouse(:, term, m), '--','Color',coloruse(term),'linewidth',0.5);
    end
    axlb(term) = plot(bin_centers, frac_sig(:,term),'Color',coloruse(term), 'linewidth',2);
end
xlabel('Time to choice (s)'); ylabel('Fraction selective');
legend(axlb, labels)

% b) choice history * temp wont work because I ran things wrong! *
labels = {'Choice','Reward','Previous Choice','Previous Reward'};
frac_sig_hist = zeros(length(bin_edges)-1, length(labels));
for term = 1:length(labels)
    num_significant = squeeze(p_mdl_hist(:, term, :) < 0.05);
    frac_sig_hist(:,term) = sum(num_significant, 2) / length(unitIDsave);
end
frac_sig_hist_mouse = zeros(length(bin_edges)-1, length(labels), length(miceUse));
for term = 1:4
    for m = 1:length(miceUse)
        num_significant = squeeze(p_mdl_hist(:, term, mouseIDsave==m) < 0.05);
        frac_sig_hist_mouse(:, term, m) = sum(num_significant, 2) / sum(mouseIDsave==m);
    end
end
% plot
figure; hold on; coloruse = 'brcm'; axlb = [];
for term = 1:length(labels)
    for m = 1:length(miceUse)
        plot(bin_centers, frac_sig_hist_mouse(:, term, m), '--','Color',coloruse(:,term),'linewidth',0.5);
    end
    axlb(term) = plot(bin_centers, frac_sig_hist(:,term),'Color',coloruse(term), 'linewidth',2);
end
xlabel('Time to choice (s)'); ylabel('Fraction selective');
legend(axlb, labels)

% c) value <- this should have longer legs than history
frac_sig_Q = zeros(length(bin_edges)-1, 1);
frac_sig_Q_mouse = zeros(length(bin_edges)-1, length(miceUse));
for m = 1:length(miceUse)
    num_significant = squeeze(p_mdl_Q(:,1,mouseIDsave==m) < 0.05);
    frac_sig_Q_mouse(:, m) = sum(num_significant, 2) / sum(mouseIDsave==m);
end
num_significant = squeeze(p_mdl_Q(:, 1, :) < 0.05);
frac_sig_Q(:, 1) = sum(num_significant, 2) / length(unitIDsave);

% Lets plot everything?
% plot
figure; hold on; coloruse = 'brg'; axlb = [];
for m = 1:length(miceUse)
    plot(bin_centers, frac_sig_mouse(:, 1, m), '--','Color',coloruse(1),'linewidth',0.5); % choice
    plot(bin_centers, frac_sig_mouse(:, 2, m), '--','Color',coloruse(2),'linewidth',0.5); % reward
    plot(bin_centers, frac_sig_Q_mouse(:, m), '--', 'Color',coloruse(3),'linewidth',0.5); % value
end
ax1 = plot(bin_centers, frac_sig(:,1),'Color',coloruse(1), 'linewidth',2);
ax2 = plot(bin_centers, frac_sig(:,2),'Color',coloruse(2), 'linewidth',2);
ax3 = plot(bin_centers, frac_sig_Q(:,1),'Color',coloruse(3), 'linewidth',2);

xlabel('Time to choice (s)'); ylabel('Fraction selective');
legend([ax1,ax2,ax3], {'Choice','Reward','Value'})

% d) hist 2 back
labels = {'Choice','Reward','Choice (t-1)','Reward (t-1)','Choice (t-2)','Reward (t-2)'};
frac_sig_hist = zeros(length(bin_edges)-1, length(labels));
for term = 1:length(labels)
    num_significant = squeeze(p_mdl_hist_2(:, term, :) < 0.05);
    frac_sig_hist(:,term) = sum(num_significant, 2) / length(unitIDsave);
end
frac_sig_hist_mouse = zeros(length(bin_edges)-1, length(labels), length(miceUse));
for term = 1:6
    for m = 1:length(miceUse)
        num_significant = squeeze(p_mdl_hist_2(:, term, mouseIDsave==m) < 0.05);
        frac_sig_hist_mouse(:, term, m) = sum(num_significant, 2) / sum(mouseIDsave==m);
    end
end
% plot
figure; hold on; coloruse = 'brcm'; axlb = [];
coloruse = [0,1,0; 1,0,0; 0, 1,1; 1,0,1; 0,0.5,0.5; 0.5,0,0.5];
for term = 1:length(labels)
    for m = 1:length(miceUse)
        plot(bin_centers, frac_sig_hist_mouse(:, term, m), '--','Color',coloruse(term,:),'linewidth',0.5);
    end
    axlb(term) = plot(bin_centers, frac_sig_hist(:,term),'Color',coloruse(term,:), 'linewidth',2);
end
xlabel('Time to choice (s)'); ylabel('Fraction selective');
legend(axlb, labels)

%% Encoding strenght plots?
% - histogram of encoding weights, comparisong social vs. nonsocial

% Change to use the balanced class version?
% average of weights across the entire time period?
nTerms = size(p_mdl_social_unbalanced,2);
weight_sig = zeros(size(p_mdl_social_unbalanced,3), size(p_mdl_social_unbalanced, 2));
for term = 1:nTerms
    num_significant = squeeze(p_mdl_social_unbalanced(:, term, :) < 0.05);    
    b = squeeze(b_mdl_social_unbalanced(:,term,:));
    % set non significant terms to 0
    b(num_significant==0) = nan;
    weight_sig(:, term) = mean(b,1, 'omitnan');
end
edges = -5:0.5:5;

figure; 
subplot(1,3,1); hold on;
histogram(weight_sig(:,1),edges,'Normalization','Probability'); 
histogram(weight_sig(:,3),edges,'Normalization','Probability');
legend('Self','Other'); title('Choice'); ylabel('Fraction');
subplot(1,3,2); hold on;
histogram(weight_sig(:,2),edges,'Normalization','Probability'); 
histogram(weight_sig(:,4),edges,'Normalization','Probability');
legend('Self','Other'); title('Reward'); ylabel('Fraction');
xlabel('Average encoding weights on significant fits');

% add value?
weight_val = zeros(size(b_mdl_Q_social,3), size(b_mdl_Q_social, 2));
for term = 1:2
    num_significant = squeeze(p_mdl_Q_social(:, term, :) < 0.05);
    b = squeeze(b_mdl_Q_social(:, term, :));
    b(num_significant==0) = nan;
    weight_val(:,term) = mean(b, 1, 'omitnan');
end
subplot(1,3,3); hold on;
histogram(weight_val(:,1),edges,'Normalization','probability');
histogram(weight_val(:,2),edges,'Normalization','probability');
legend('Self','Other'); title('Value'); ylabel('Fraction');

% CDF plots
figure; 
ax1=subplot(1,3,1); hold on;
cdfplot(weight_sig(:,1)); cdfplot(weight_sig(:,3));
legend('Self','Other'); xlabel('\beta'); ylabel('Cumulative fraction');
title('Choice');
ax2=subplot(1,3,2); hold on;
cdfplot(weight_sig(:,2)); cdfplot(weight_sig(:,4));
legend('Self','Other'); xlabel('\beta'); ylabel('Cumulative fraction');
title('Reward');
ax3=subplot(1,3,3); hold on;
cdfplot(weight_val(:,1)); cdfplot(weight_val(:,2));
legend('Self','Other'); xlabel('\beta'); ylabel('Cumulative fraction');
title('Value')
linkaxes([ax1,ax2,ax3],'xy');

%% Weight comparison for joint model

% for a given time window?
twin = 6:26; % uh?

twin = 1:length(bin_centers);

nTerms = size(p_mdl_social_joint,2);
weight_sig = zeros(size(p_mdl_social_joint,3), size(p_mdl_social_joint, 2));
for term = 1:nTerms
    num_significant = squeeze(p_mdl_social_joint(twin, term, :) < 0.05);    
    b = squeeze(b_mdl_social_joint(twin,term,:));
    % set non significant terms to 0
    b(num_significant==0) = nan;
    weight_sig(:, term) = mean(b,1, 'omitnan');
end

edges = -3:0.25:3;

% Histograms
figure; 
subplot(1,3,1); hold on;
histogram(weight_sig(:,1),edges,'Normalization','Probability','EdgeColor','none','FaceAlpha',0.5); 
histogram(weight_sig(:,4),edges,'Normalization','Probability','EdgeColor','none','FaceAlpha',0.5);
legend('Self','Other'); title('Choice'); ylabel('Fraction');
subplot(1,3,2); hold on;
histogram(weight_sig(:,2),edges,'Normalization','Probability','EdgeColor','none','FaceAlpha',0.5); 
histogram(weight_sig(:,5),edges,'Normalization','Probability','EdgeColor','none','FaceAlpha',0.5);
legend('Self','Other'); title('Reward'); ylabel('Fraction');
xlabel('Average encoding weights on significant fits');
subplot(1,3,3); hold on;
histogram(weight_sig(:,3),edges,'Normalization','Probability','EdgeColor','none','FaceAlpha',0.5); 
histogram(weight_sig(:,6),edges,'Normalization','Probability','EdgeColor','none','FaceAlpha',0.5);
legend('Self','Other'); title('Value'); ylabel('Fraction');

[~,p_choice] = kstest2(weight_sig(:,1), weight_sig(:,4));
[~,p_reward] = kstest2(weight_sig(:,2), weight_sig(:,5));
[~,p_value] = kstest2(weight_sig(:,3), weight_sig(:,6));

% CDF
figure; 
ax1=subplot(1,3,1); hold on;
cdfplot(weight_sig(:,1)); cdfplot(weight_sig(:,4));
legend('Self','Other'); xlabel('\beta'); ylabel('Cumulative fraction');
title(['Choice, p_ks = ' num2str(p_choice)]);
ax2=subplot(1,3,2); hold on;
cdfplot(weight_sig(:,2)); cdfplot(weight_sig(:,5));
legend('Self','Other'); xlabel('\beta'); ylabel('Cumulative fraction');
title(['Reward, p_ks = ' num2str(p_reward)]);
ax3=subplot(1,3,3); hold on;
cdfplot(weight_sig(:,3)); cdfplot(weight_sig(:,6));
legend('Self','Other'); xlabel('\beta'); ylabel('Cumulative fraction');
title(['Value, p_ks = ' num2str(p_value)])
linkaxes([ax1,ax2,ax3],'xy');

xlim([-5,5])



%% Social plots

p_mdl_social_unbalanced; % chocie, reward, chocie, reward

nTerms = size(p_mdl_social_unbalanced,2);
frac_sig = zeros(length(bin_edges)-1, nTerms);
frac_sig_mouse = zeros(length(bin_edges)-1, nTerms, length(miceUse));
for term = 1:nTerms
    num_significant = squeeze(p_mdl_social_unbalanced(:, term, :) < 0.05);
    frac_sig(:,term) = sum(num_significant, 2) / length(unitIDsave);
    for m = 1:length(miceUse)
        num_significant = squeeze(p_mdl_social_unbalanced(:, term, mouseIDsave==m) < 0.05);
        frac_sig_mouse(:, term, m) = sum(num_significant, 2) / sum(mouseIDsave==m);
    end
end
% plot
figure; hold on; 
coloruse = 'brg';
coloruse = [0.5,0.5,1,0.5; 1,0.5,0.5,0.5; 0,0,0.7,1; 0.7,0,0,1];
axlb = [];
for term = 1:nTerms
    for m = 1:length(miceUse)
        plot(bin_centers, frac_sig_mouse(:, term, m), '--','Color',coloruse(term,:),'linewidth',0.5); 
    end
    axlb(term) = plot(bin_centers, frac_sig(:, term), 'Color',coloruse(term,:), 'linewidth',2);
end
legend(axlb,{'Choice','Reward','Social choice','Social Reward'});
xlabel('Time to choice (s)'); ylabel('Fraction significant neurons');

% Include value?
% c) value <- this should have longer legs than history
nTerms = 2;
frac_sig_Q = zeros(length(bin_edges)-1, nTerms);
frac_sig_Q_mouse = zeros(length(bin_edges)-1, nTerms, length(miceUse));
for term = 1:nTerms
for m = 1:length(miceUse)
    num_significant = squeeze(p_mdl_Q_social(:,term,mouseIDsave==m) < 0.05);
    frac_sig_Q_mouse(:, term, m) = sum(num_significant, 2) / sum(mouseIDsave==m);
end
num_significant = squeeze(p_mdl_Q_social(:, term, :) < 0.05);
frac_sig_Q(:, term) = sum(num_significant, 2) / length(unitIDsave);
end

coloruse = [0.5, 1, 0.5, 0.5; 0, 0.7, 0,1];
for term = 1:nTerms
    for m = 1:length(miceUse)
        plot(bin_centers, frac_sig_Q_mouse(:, term, m), '--', 'Color', coloruse(term,:), 'linewidth',0.5 )
    end
    axlb(end+1) = plot(bin_centers, frac_sig_Q(:, term), 'Color', coloruse(term,:), 'linewidth', 2);
end

legend(axlb,{'Choice','Reward','Social choice','Social Reward','Value','Social Value'});
xlabel('Time to choice (s)'); ylabel('Fraction significant neurons');


%% Plot with balanced classes?
% choice, reward, choice, reward
frac_sig_balanced_choice = [];
num_significant = squeeze(p_mdl_social_balanced_choice(:, 1, :, :) < 0.05);
frac_sig_balanced_choice(:,1) = sum(num_significant, [2,3]) ./ (length(unitIDsave) * nShuff);
num_significant = squeeze(p_mdl_social_balanced_choice(:, 3, :, :) < 0.05);
frac_sig_balanced_choice(:,2) = sum(num_significant, [2,3]) ./ (length(unitIDsave) * nShuff);

frac_sig_balanced_choice_mouse = [];
for m = 1:length(miceUse)
    num_significant = squeeze(p_mdl_social_balanced_choice(:, 1, mouseIDsave==m, :) < 0.05);
    frac_sig_balanced_choice_mouse(:, 1, m) = sum(num_significant, [2,3]) / (sum(mouseIDsave==m) * nShuff );
    num_significant = squeeze(p_mdl_social_balanced_choice(:, 3, mouseIDsave==m, :) < 0.05);
    frac_sig_balanced_choice_mouse(:, 2, m) = sum(num_significant, [2,3]) / (sum(mouseIDsave==m) * nShuff );
end

frac_sig_balanced_reward = [];
num_significant = squeeze(p_mdl_social_balanced_reward(:, 2, :, :) < 0.05);
frac_sig_balanced_reward(:,1) = sum(num_significant, [2,3]) ./ (length(unitIDsave) * nShuff);
num_significant = squeeze(p_mdl_social_balanced_reward(:, 4, :, :) < 0.05);
frac_sig_balanced_reward(:,2) = sum(num_significant, [2,3]) ./ (length(unitIDsave) * nShuff);

frac_sig_balanced_reward_mouse = [];
for m = 1:length(miceUse)
    num_significant = squeeze(p_mdl_social_balanced_reward(:, 2, mouseIDsave==m, :) < 0.05);
    frac_sig_balanced_reward_mouse(:, 1, m) = sum(num_significant, [2,3]) / (sum(mouseIDsave==m) * nShuff );
    num_significant = squeeze(p_mdl_social_balanced_reward(:, 4, mouseIDsave==m, :) < 0.05);
    frac_sig_balanced_reward_mouse(:, 2, m) = sum(num_significant, [2,3]) / (sum(mouseIDsave==m) * nShuff );
end

figure; hold on;

% plot the previous social only?
coloruse = 'brg';
coloruse = [0.5,0.5,1; 1,0.5,0.5; 0,0,0.7; 0.7,0,0];
axlb = [];
for term = 3:4
    for m = 1:length(miceUse)
        plot(bin_centers, frac_sig_mouse(:, term, m), '--','Color',coloruse(term,:),'linewidth',0.5); 
    end
    axlb(end+1) = plot(bin_centers, frac_sig(:, term), 'Color',coloruse(term,:), 'linewidth',2);
end
xlabel('Time to choice (s)'); ylabel('Fraction significant neurons');


for m = 1:length(miceUse)
    plot(bin_centers, frac_sig_balanced_choice_mouse(:, 2, m), ':','Color',[0,0,1],'linewidth',0.5)
end
axlb(end+1) = plot(bin_centers, frac_sig_balanced_choice(:, 2),'Color',[0,0,1],'linewidth',2);

for m = 1:length(miceUse)
    plot(bin_centers, frac_sig_balanced_reward_mouse(:, 2, m), ':','Color',[1,0,0],'linewidth',0.5)
end
axlb(end+1) = plot(bin_centers, frac_sig_balanced_reward(:, 2),'Color',[1,0,0],'linewidth',2);

legend(axlb,{'Social choice','Social Reward','Social choice - balanced','Social Reward - balanced'});


%% Plot with social shuffle?


p_mdl_social_shuffled;

% calculate fraction significant for shuffled
frac_sig_shuffled = [];
frac_sig_shuffled_mouse = [];
for term = 1:4
for m = 1:length(miceUse)
    num_significant = squeeze(p_mdl_social_shuffled(:, term, mouseIDsave==m, :) < 0.05);
    frac_sig_shuffled_mouse(:, term, m) = sum(num_significant, [2,3]) / (sum(mouseIDsave==m) * nShuff );
    
end
num_significant = squeeze(p_mdl_social_shuffled(:, term, :, :) < 0.05);
frac_sig_shuffled(:,term) = sum(num_significant, [2,3]) ./ (length(unitIDsave) * nShuff);
end
% shoudl use control here 
p_mdl_social_shuffled_control;
frac_sig_shuffled_control_mouse = [];
frac_sig_shuffled_control = [];
for term = 1:4
for m = 1:length(miceUse)
    num_significant = squeeze(p_mdl_social_shuffled_control(:, term, mouseIDsave==m) < 0.05);
    frac_sig_shuffled_control_mouse(:, term, m) = sum(num_significant, [2,3]) / (sum(mouseIDsave==m) );
    
end
num_significant = squeeze(p_mdl_social_shuffled_control(:, term, :) < 0.05);
frac_sig_shuffled_control(:,term) = sum(num_significant, [2,3]) ./ (length(unitIDsave));
end

figure; hold on;
% plot the previous social only?
coloruse = 'brg';
coloruse = [0.5,0.5,1; 1,0.5,0.5; 0,0,0.7; 0.7,0,0];
axlb = [];
for term = 3:4
    for m = 1:length(miceUse)
        %plot(bin_centers, frac_sig_mouse(:, term, m), '--','Color',coloruse(term,:),'linewidth',0.5); 
        plot(bin_centers, frac_sig_shuffled_control_mouse(:, term, m), '--','Color',coloruse(term,:),'linewidth',0.5); 
    end
    %axlb(end+1) = plot(bin_centers, frac_sig(:, term), 'Color',coloruse(term,:), 'linewidth',2);
    axlb(end+1) = plot(bin_centers, frac_sig_shuffled_control(:, term), 'Color',coloruse(term,:), 'linewidth',2);
end
xlabel('Time to choice (s)'); ylabel('Fraction significant neurons');

% plot the shuffle control
coloruse = [0.5,0.5,1; 1,0.5,0.5; 0,0,1; 1,0,0];
for term = 3:4
    for m = 1:length(miceUse)
        plot(bin_centers, frac_sig_shuffled_mouse(:, term, m), ':','Color',coloruse(term,:), 'linewidth',0.5 )
    end
    axlb(end+1) = plot(bin_centers, frac_sig_shuffled(:, term), 'Color', coloruse(term,:), 'linewidth',2);
end

legend(axlb,{'Social choice','Social Reward','Social choice - shuffled','Social Reward - shuffled'});


%% Value shuffle

p_mdl_Q_social_shuffled;

% calculate fraction significant for shuffled
frac_sig_Q_shuffled = [];
frac_sig_Q_shuffled_mouse = [];
for term = 1:2
for m = 1:length(miceUse)
    num_significant = squeeze(p_mdl_Q_social_shuffled(:, term, mouseIDsave==m, :) < 0.05);
    frac_sig_Q_shuffled_mouse(:, term, m) = sum(num_significant, [2,3]) / (sum(mouseIDsave==m) * nShuff );
    
end
num_significant = squeeze(p_mdl_Q_social_shuffled(:, term, :, :) < 0.05);
frac_sig_Q_shuffled(:,term) = sum(num_significant, [2,3]) ./ (length(unitIDsave) * nShuff);
end
% shoudl use control here 
p_mdl_Q_social_shuffled_control;
frac_sig_Q_shuffled_control_mouse = [];
frac_sig_Q_shuffled_control = [];
for term = 1:2
for m = 1:length(miceUse)
    num_significant = squeeze(p_mdl_Q_social_shuffled_control(:, term, mouseIDsave==m) < 0.05);
    frac_sig_Q_shuffled_control_mouse(:, term, m) = sum(num_significant, [2,3]) / (sum(mouseIDsave==m) );
    
end
num_significant = squeeze(p_mdl_Q_social_shuffled_control(:, term, :) < 0.05);
frac_sig_Q_shuffled_control(:,term) = sum(num_significant, [2,3]) ./ (length(unitIDsave));
end

figure; hold on;
% plot the previous social only?
coloruse = 'brg';
coloruse = [0.5,1,0.5;0,0.7,0];
axlb = [];
for term = 2
    for m = 1:length(miceUse)
        plot(bin_centers, frac_sig_Q_shuffled_control_mouse(:, term, m), '--','Color',coloruse(term,:),'linewidth',0.5); 
    end
    axlb(end+1) = plot(bin_centers, frac_sig_Q_shuffled_control(:, term), 'Color',coloruse(term,:), 'linewidth',2);
end
xlabel('Time to choice (s)'); ylabel('Fraction significant neurons');

% plot the shuffle control
coloruse = [0.5,1,0.5,; 0,1,0];
for term = 2
    for m = 1:length(miceUse)
        plot(bin_centers, frac_sig_Q_shuffled_mouse(:, term, m), ':','Color',coloruse(term,:), 'linewidth',0.5 )
    end
    axlb(end+1) = plot(bin_centers, frac_sig_Q_shuffled(:, term), 'Color', coloruse(term,:), 'linewidth',2);
end

legend(axlb,{'Value','Value shuffled'});

%% Social plot with joint encoding


nTerms = size(p_mdl_social_joint,2);
frac_sig = zeros(length(bin_edges)-1, nTerms);
frac_sig_mouse = zeros(length(bin_edges)-1, nTerms, length(miceUse));
for term = 1:nTerms
    num_significant = squeeze(p_mdl_social_joint(:, term, :) < 0.05);
    frac_sig(:,term) = sum(num_significant, 2) / length(unitIDsave);
    for m = 1:length(miceUse)
        num_significant = squeeze(p_mdl_social_joint(:, term, mouseIDsave==m) < 0.05);
        frac_sig_mouse(:, term, m) = sum(num_significant, 2) / sum(mouseIDsave==m);
    end
end
% plot
figure; hold on; 
coloruse = 'brg';
coloruse = [0.5,0.5,1,0.5; 1,0.5,0.5,0.5;  0.5, 1, 0.5, 0.5 ;...
    0,0,0.7,1; 0.7,0,0,1; 0,0.7,0,1];
axlb = [];
for term = 1:nTerms
    for m = 1:length(miceUse)
        plot(bin_centers, frac_sig_mouse(:, term, m), '--','Color',coloruse(term,:),'linewidth',0.5); 
    end
    axlb(term) = plot(bin_centers, frac_sig(:, term), 'Color',coloruse(term,:), 'linewidth',2);
end
legend(axlb,{'Choice','Reward','Value','Social choice','Social Reward','Social Value'});
xlabel('Time to choice (s)'); ylabel('Fraction significant neurons');

% non social joint
nTerms = size(p_mdl_joint,2);
frac_sig = zeros(length(bin_edges)-1, nTerms);
frac_sig_mouse = zeros(length(bin_edges)-1, nTerms, length(miceUse));
for term = 1:nTerms
    num_significant = squeeze(p_mdl_joint(:, term, :) < 0.05);
    frac_sig(:,term) = sum(num_significant, 2) / length(unitIDsave);
    for m = 1:length(miceUse)
        num_significant = squeeze(p_mdl_joint(:, term, mouseIDsave==m) < 0.05);
        frac_sig_mouse(:, term, m) = sum(num_significant, 2) / sum(mouseIDsave==m);
    end
end
% plot
figure; hold on; 
coloruse = 'brg';
coloruse = [0,0,1; 1,0,0,; 0,1,0];
axlb = [];
for term = 1:nTerms
    for m = 1:length(miceUse)
        plot(bin_centers, frac_sig_mouse(:, term, m), '--','Color',coloruse(term,:),'linewidth',0.5); 
    end
    axlb(term) = plot(bin_centers, frac_sig(:, term), 'Color',coloruse(term,:), 'linewidth',2);
end
legend(axlb,{'Choice','Reward','Value'});
xlabel('Time to choice (s)'); ylabel('Fraction significant neurons');


%% smaller checks

%% Does social encoding decrease with the training time? Or is this just a function of not having enough samples?
% Yes I believe this is a function of decreasing sample size...
% - need to balance classes and sample size.
% - issue is mice get better and do fewer opposite trial runs

frac_sig_balanced_choice_mouse_date = [];
frac_sig_unbalanced_choice_mouse_date = [];
for m = 1:length(miceUse)
    for d = 1:length(dateUse{m})
        iduse = find(mouseIDsave==m & dateIDsave==d);
        num_significant = squeeze(p_mdl_social_balanced_choice(:, 1, iduse, :) < 0.05);
        frac_sig_balanced_choice_mouse_date(:, 1, m,d) = sum(num_significant, [2,3]) / (length(iduse) * nShuff );
        num_significant = squeeze(p_mdl_social_balanced_choice(:, 3, iduse, :) < 0.05);
        frac_sig_balanced_choice_mouse_date(:, 2, m,d) = sum(num_significant, [2,3]) / (length(iduse) * nShuff );

        num_significant = squeeze(p_mdl_social_unbalanced(:, 1, iduse) < 0.05);
        frac_sig_unbalanced_choice_mouse_date(:, 1, m, d) = sum(num_significant, 2) / length(iduse);
        num_significant = squeeze(p_mdl_social_unbalanced(:, 3, iduse) < 0.05);
        frac_sig_unbalanced_choice_mouse_date(:, 2, m, d) = sum(num_significant, 2) / length(iduse);
    end
end

m = 3;
figure; hold on;
for d = 1:3
    plot(bin_centers, frac_sig_balanced_choice_mouse_date(:, 1, m, d))
    %plot(bin_centers, frac_sig_unbalanced_choice_mouse_date(:, 2, m, d))
end
legend('Day 1','Day 2', 'Day 3')



%% Plot balanced classes for non social

figure; hold on;

% plot the previous social only?
coloruse = 'brg';
coloruse = [0.5,0.5,1; 1,0.5,0.5; 0,0,1; 1,0,0];
axlb = [];
for term = 1:2
    for m = 1:length(miceUse)
        plot(bin_centers, frac_sig_mouse(:, term, m), '--','Color',coloruse(term,:),'linewidth',0.5); 
    end
    axlb(end+1) = plot(bin_centers, frac_sig(:, term), 'Color',coloruse(term,:), 'linewidth',2);
end
xlabel('Time to choice (s)'); ylabel('Fraction significant neurons');


for m = 1:length(miceUse)
    plot(bin_centers, frac_sig_balanced_choice_mouse(:, 1, m), ':','Color',[0,0,0.6],'linewidth',0.5)
end
axlb(end+1) = plot(bin_centers, frac_sig_balanced_choice(:, 1),'Color',[0,0,0.6],'linewidth',2);

for m = 1:length(miceUse)
    plot(bin_centers, frac_sig_balanced_reward_mouse(:, 1, m), ':','Color',[0.6,0,0],'linewidth',0.5)
end
axlb(end+1) = plot(bin_centers, frac_sig_balanced_reward(:, 1),'Color',[0.6,0,0],'linewidth',2);

legend(axlb,{' choice',' Reward',' choice - balanced',' Reward - balanced'});

%% Plot interaction terms!
% e.g., choice1 x choice2

interaction_labels_social = {'Choice','Reward','Interaction','Social choice','Social reward',...
    'Social interaction','Joint choice','Joint reward','Joint everything!'};

nTerms = size(p_mdl_interaction_social,2);
frac_sig = zeros(length(bin_edges)-1, nTerms);
frac_sig_mouse = zeros(length(bin_edges)-1, nTerms, length(miceUse));
for term = 1:nTerms
    num_significant = squeeze(p_mdl_interaction_social(:, term, :) < 0.05);
    frac_sig(:,term) = sum(num_significant, 2) / length(unitIDsave);
    for m = 1:length(miceUse)
        num_significant = squeeze(p_mdl_interaction_social(:, term, mouseIDsave==m) < 0.05);
        frac_sig_mouse(:, term, m) = sum(num_significant, 2) / sum(mouseIDsave==m);
    end
end
% plot
figure; hold on; 
coloruse = linspecer(nTerms);
axlb = [];
for term = 1:nTerms
    for m = 1:length(miceUse)
        plot(bin_centers, frac_sig_mouse(:, term, m), '--','Color',coloruse(term,:),'linewidth',0.5); 
    end
    axlb(term) = plot(bin_centers, frac_sig(:, term), 'Color',coloruse(term,:), 'linewidth',2);
end
legend(axlb,interaction_labels_social);
xlabel('Time to choice (s)'); ylabel('Fraction significant neurons');

% non social joint
nTerms = size(p_mdl_interaction,2);
frac_sig = zeros(length(bin_edges)-1, nTerms);
frac_sig_mouse = zeros(length(bin_edges)-1, nTerms, length(miceUse));
for term = 1:nTerms
    num_significant = squeeze(p_mdl_interaction(:, term, :) < 0.05);
    frac_sig(:,term) = sum(num_significant, 2) / length(unitIDsave);
    for m = 1:length(miceUse)
        num_significant = squeeze(p_mdl_interaction(:, term, mouseIDsave==m) < 0.05);
        frac_sig_mouse(:, term, m) = sum(num_significant, 2) / sum(mouseIDsave==m);
    end
end
% plot
figure; hold on; 
coloruse = 'brg';
coloruse = [0,0,1; 1,0,0,; 0,1,0];
axlb = [];
for term = 1:nTerms
    for m = 1:length(miceUse)
        plot(bin_centers, frac_sig_mouse(:, term, m), '--','Color',coloruse(term,:),'linewidth',0.5); 
    end
    axlb(term) = plot(bin_centers, frac_sig(:, term), 'Color',coloruse(term,:), 'linewidth',2);
end
legend(axlb,{'Choice','Reward','Choice x Reward'});
xlabel('Time to choice (s)'); ylabel('Fraction significant neurons');


%% 2) 

%% Mixed selectivity
% Non-social

time_window_crop = [-2, 2];
twin = 6:26; % -2 to 2 seconds, emperically

model_labels = {'C','R'};

minBins = 3; % bin_size = 0.2 seconds, so >0.6 seconds in a row its significant
[~, num_models, num_neurons] = size(p_mdl);
sig_models = zeros(num_models, num_neurons);

sig_models_labels = cell(num_neurons,1);

for n = 1:num_neurons
    for m = 1:num_models
        p = p_mdl(twin, m, n) < 0.05;
        sig = conv(double(p), ones(minBins, 1), 'valid')==minBins;
        sig_models(m, n) = any(sig);
    end
    tmp = model_labels(sig_models(:,n)==1);
    sig_models_labels{n} = [tmp{:}];
end
sig_models_labels(cellfun(@isempty, sig_models_labels)) = {'None'};

figure;
histogram(categorical(sig_models_labels))
ylabel('Count');
hold on; plot(xlim, [1,1] * num_neurons * 0.05, 'k--');
title('Tuning over time bins')

% for social
if contains(mouse_name, '-')
    model_labels = {'C','R','c','r'};

    [~, num_models, num_neurons] = size(p_mdl_social_unbalanced);
    sig_models = zeros(num_models, num_neurons);
    
    sig_models_labels = cell(num_neurons,1);
    
    for n = 1:num_neurons
        for m = 1:num_models
            p = p_mdl_social_unbalanced(twin, m, n) < 0.05;
            sig = conv(double(p), ones(minBins, 1), 'valid')==minBins;
            sig_models(m, n) = any(sig);
        end
        tmp = model_labels(sig_models(:,n)==1);
        sig_models_labels{n} = [tmp{:}];
    end
    sig_models_labels(cellfun(@isempty, sig_models_labels)) = {'None'};
    
    figure;
    histogram(categorical(sig_models_labels))
    ylabel('Count');
    hold on; plot(xlim, [1,1] * num_neurons * 0.05, 'k--');

end

%% plot social categories
% only looking within a category
p_mdl_social = cat(2,p_mdl_social_unbalanced(:, [1,2], :) , p_mdl_Q_social(:, 1, :), ...
    p_mdl_social_unbalanced(:, [3,4], :) , p_mdl_Q_social(:, 2, :));

p_mdl_social =  p_mdl_social_joint;

proportions_saved = []; counts_save = [];

ver_plot_use = 2;
ax1 = [];

figure;
terms_use_loop = {[1,4],[2,5],[3,6]};
title_label = {'Choice','Reward','Value'};
for j = 1:3
terms_use = terms_use_loop{j};
if contains(mouse_name, '-')
    model_labels = {'C','R','V','c','r','v'};


    model_labels = model_labels(terms_use);

    [~, num_models, num_neurons] = size(p_mdl_social);
    sig_models = zeros(num_models, num_neurons);
    
    sig_models_labels = cell(num_neurons,1);
    
    for n = 1:num_neurons
        for m = 1:num_models
            p = p_mdl_social(twin, m, n) < 0.05;
            %if m==4; p = p(1:16); end % temp crops out after choice?
            sig = conv(double(p), ones(minBins, 1), 'valid')==minBins;
            sig_models(m, n) = any(sig);
        end
        tmp = model_labels(sig_models(terms_use,n)==1);
        sig_models_labels{n} = [tmp{:}];
    end
    sig_models_labels(cellfun(@isempty, sig_models_labels)) = {'None'};
    
    %figure;
    if ver_plot_use==1
    subplot(1,3,j);
    histogram(categorical(sig_models_labels))
    ylabel('Count');
    hold on; plot(xlim, [1,1] * num_neurons * 0.05, 'k--');
    %hold on; plot(xlim, [1,1]*0.05, 'k--')
    title(title_label{j})


    elseif ver_plot_use==2

    % different plot without none and normalized
    ax1(j) = subplot(1,3,j); hold on;
    % Get unique categories and their counts
    [cats, ~, idx] = unique(sig_models_labels);
    counts = accumarray(idx, 1);
    norm_counts = counts / length(sig_models_labels);
    
    % Remove 'None' from categories and normalized counts
    keep = ~strcmp(cats, 'None');  % logical index of categories to keep
    
    % Final data to plot
    plot_cats = cats(keep);
    plot_vals = norm_counts(keep);
    
    % Bar plot
    bar(categorical(plot_cats), plot_vals);
    ylabel('Normalized Frequency');
    title(title_label{j})
    hold on; plot(xlim, [1,1] * 0.05, 'k--');

    % Statistics:
    %p1 = 

    % Save:
    proportions_save(:, j) = plot_vals;
    counts_save(:,j) = counts(keep);


    end
end
end
linkaxes(ax1, 'y')

%% Same as above but as a venn diagram

% get scale factors
totals = sum(counts_save,1);


term_labels = {'Choice','Reward','Value'};
% uh
figure; 
for term_use = 1:3
Nonsocial = counts_save(1,term_use);
Social = counts_save(3, term_use);
Joint = counts_save(2, term_use);

% Set sizes
A = Nonsocial + Joint;   % Size of Set A
B = Social + Joint;   % Size of Set B
AB = Joint;  % Intersection size (A ∩ B)

scale = sqrt(totals(term_use) / max(totals));
%scale = (totals(term_use) / num_neurons);

subplot(1,3,term_use);
plotVennDiagramProportional(A, B, AB, scale);
title(term_labels{term_use});
end

%% Use joint model for the above?


