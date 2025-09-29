%% Generate basic learning plots to determine expert level

addpath('../')


miceNonsocial = {'KM35','KM36','KM37','KM38'};
miceSocial = {'KM35-36','KM37-38'};
% % 
% miceNonsocial = {'KM33','KM34'};
% miceSocial = {'KM33-34'};
% % 
%   miceNonsocial = {'KM39','KM40','KM41','KM42'};
%   miceSocial = {'KM39-40','KM41-42'};
% 
% miceNonsocial = {'KM43','KM44'};
% miceNonsocial = {'KM45','KM46','KM47','KM48'};
% 
% miceSocial = {};

miceNonsocial = {'KM10','KM12','KM11','KM13','KM20','KM23','KM21','KM22','KM24','KM26','KM25','KM27','KM28','KM29','KM30','KM31'};
miceSocial = {'KM10-12','KM11-13','KM20-23','KM21-22','KM24-26','KM25-27','KM28-29','KM30-31'};

miceNonsocial = {'KM49','KM50','KM51','KM52','KM53','KM54','KM55','KM56'};
miceSocial = {};

miceNonsocial = {'KM49','KM50'};
miceSocial = {'KM49-50'};

use_social_minus = 0;

do_pad = 1; % ?

data_path = 'Z:\HarveyLab\Tier1\Kevin\Videos\';

username = 'KM';

%%
% data_path = 'Z:\HarveyLab\Tier1\Rhyanne';
% % % miceNonsocial = {'RF1','RF2','RF3','RF4','RF5','RF6','RF7'};
% % % groupID = [0,0,1,0,0,0,1];
% % % miceSocial = {};
% % 
% miceNonsocial = {'RF3','RF4','RF5','RF6','RF7'};
% miceSocial = {'RF4-6','RF5-7'};
% groupID = [1,0,0,0,1];
% % 
% miceNonsocial = {'RF8','RF9','RF10','RF11','RF12','RF13','RF14'};
% miceSocial = {'RF8-10','RF11-12','RF13-14'};
% groupID = [1,0,0, 0,1,0,0];
% % 
% % miceNonsocial = {'RF8','RF9','RF10'};
% % miceSocial = {};
% % groupID = [1,0,0];
% % % miceNonsocial = {'RF11','RF12','RF13','RF14'};
% % % miceSocial = {};
% % % groupID = [ 0,1,0,0];
% 
% 
% username = 'RF';

%% Non social
if ~exist('groupID','var'); groupID = []; end

%prob_use_mice = 80 * ones(1,length(miceNonsocial));
if ~isequal(size(groupID), size(miceNonsocial))
    groupID = zeros(size(miceNonsocial));
end

doPad = 1; % pad instead of crop

accuracy_per_sess = {};
trials_per_sess = {}; % or trials per min? I think this is better...
trial_time_sess = {}; % this might be using stem time, not door open time...
sesstime_all = {};
% start time not trial time

% new, look at trial times and how its contingent on accuracy?
acc_given_ITI = {}; % <- this is more meaningful
acc_given_choicetime = {}; % <- this is meaningless, since I have random delays!!!
acc_given_ITI_hit = {};
acc_given_ITI_miss = {};

nMice = length(miceNonsocial);

ITI_all = cell(nMice,1); acc_prev = cell(nMice,1);
protocol_save = {};
for m = 1:length(miceNonsocial)

    mouse_name = miceNonsocial{m};
    datapath = fullfile(data_path, mouse_name, 'mouseBEHstruct.mat');
    load(datapath);
    
    for s = 1:length(info)%(length(info)-sess_use) : length(info)%1:length(info)
        LTProb = info(s).LTProb;
        RTProb = info(s).RTProb;
        choice = info(s).choice;
        reward = info(s).reward;
        choice_time = info(s).choice_time;
        start_time = info(s).start_time;
        badid = cellfun(@isempty, choice) | cellfun(@isempty, reward) ...
            | cellfun(@isempty, choice_time);
        LTProb(badid) = []; choice(badid) = []; reward(badid) = []; RTProb(badid) = [];
        choice_time(badid) = []; start_time(badid) = [];
        choice = cellfun(@(v) v(1), choice);
        choice_time = cellfun(@(v) v(1), choice_time);

        %if ~any(LTProb==Prob_use(m)) || length(unique(LTProb))==1
        %if ~any(LTProb<95) || length(unique(LTProb))==1
        if length(unique([LTProb, RTProb]))==1
            accuracy_per_sess{m}(s) = nan;
            trials_per_sess{m}(s) = nan;
            trial_time_sess{m}(s) = nan;
            acc_given_ITI{m}(s,[1,2]) = nan;
            acc_given_ITI_hit{m}(s) = nan;
            acc_given_ITI_miss{m}(s) = nan;
        
        else 
            % save everything?
        accuracy = (choice>=2 & LTProb > 50) | (choice<=1 & LTProb<50);
        
        
        % accuracy over sessions
        accuracy_per_sess{m}(s) = mean(accuracy);
        
        % number of trials over sessions
        trials_per_sess{m}(s) = length(accuracy);
    
        % trial time over sessions
        trial_time = choice_time - start_time;
        trial_time_sess{m}(s) = mean(trial_time) / 1000;

        sesstime_all{m}(s) = info(s).sessionTime;

        % new, save the accuracy contingent on trial times
        start_time = info(s).start_time;
        choice_time = info(s).choice_time;
        
        choice_time(badid) = []; start_time(badid) = [];
        choice_time = cellfun(@(v) v(1), choice_time);
        
        trial_time = (choice_time - start_time) / 1000;
        %ITI = (start_time(2:end) - choice_time(1:end-1))/1000;
        % other way of calculating ITI
        ITI = (start_time(2:end) - start_time(1:end-1))/1000;
        accuracy_prev = accuracy(2:end);

        ITI_all{m} = [ITI_all{m} , ITI];
        acc_prev{m} = [acc_prev{m}, accuracy_prev];

        % save new stuff
        acc_given_ITI_hit{m}(s) = mean(ITI(accuracy_prev==1));
        acc_given_ITI_miss{m}(s) = mean(ITI(accuracy_prev==0));
        %acc_given_choicetime{m}(s,:) = [mean(trial_time(accuracy==1)), mean(trial_time(accuracy==0))];
        
        protocol_save{m}(s) = info(s).protocolNum;  
        end
    
    end
end

ITI_nonsocial = ITI_all;
acc_prev_nonsocial = acc_prev;

% rearrange
for m = 1:length(miceNonsocial)
    % accuracy_per_sess{m}(isnan(accuracy_per_sess{m})) = [];
    % trials_per_sess{m}(isnan(trials_per_sess{m})) = [];
    % trial_time_sess{m}(isnan(trial_time_sess{m})) = [];
    % acc_given_ITI_hit{m}(isnan(acc_given_ITI_hit{m})) = [];
    % acc_given_ITI_miss{m}(isnan(acc_given_ITI_miss{m})) = [];

    %*** this was 2? 
    accuracy_per_sess{m}(protocol_save{m}<1) = [];
    trials_per_sess{m}(protocol_save{m}<1) = [];
    trial_time_sess{m}(protocol_save{m}<1) = [];
    acc_given_ITI_hit{m}(protocol_save{m}<1) = [];
    acc_given_ITI_miss{m}(protocol_save{m}<1) = [];
end

accuracy_cell = accuracy_per_sess;
trials_cell = trials_per_sess;

% Pad
if doPad
    nn = max(cellfun(@length, accuracy_per_sess));
    accuracy_per_sess = cellfun(@(v) [v, nan(1, nn - length(v))], accuracy_per_sess, 'UniformOutput', false);
    trials_per_sess = cellfun(@(v) [v, nan(1, nn - length(v))], trials_per_sess, 'UniformOutput', false);
    trial_time_sess = cellfun(@(v) [v, nan(1, nn - length(v))], trial_time_sess, 'UniformOutput', false);
end

% Crop
nn = min(cellfun(@length, accuracy_per_sess));
accuracy_per_sess = cellfun(@(v) v(1:nn), accuracy_per_sess,'un',0);
accuracy_per_sess = cell2mat(accuracy_per_sess');
nn = min(cellfun(@length, trials_per_sess));
trials_per_sess = cellfun(@(v) v(1:nn), trials_per_sess,'un',0);
trials_per_sess = cell2mat(trials_per_sess');
nn = min(cellfun(@length, trial_time_sess));
trial_time_sess = cellfun(@(v) v(1:nn), trial_time_sess,'un',0);
trial_time_sess = cell2mat(trial_time_sess');
nn = min(cellfun(@length, acc_given_ITI_hit));
acc_given_ITI_hit = cellfun(@(v) v(1:nn), acc_given_ITI_hit,'un',0);
acc_given_ITI_hit = cell2mat(acc_given_ITI_hit');
acc_given_ITI_miss = cellfun(@(v) v(1:nn), acc_given_ITI_miss,'un',0);
acc_given_ITI_miss = cell2mat(acc_given_ITI_miss');



% 
accuracy_per_sess_nonsocial = accuracy_per_sess;
trials_per_sess_nonsocial = trials_per_sess;% / 60; % why divided by 60???


figure; 
ax1=subplot(1,3,1); hold on;
plot(nanmean(accuracy_per_sess,1), 'k', 'linewidth',2);
lg1 = plot(accuracy_per_sess', 'Color',[0.5,0.5,0.5,0.5],'linewidth',0.5);
lg2 = plot(accuracy_per_sess(groupID==1,:)', 'Color', [.9,.2,.2,.3],'linewidth',0.5);
ylabel('Accuracy'); xlabel('Session');

ax2=subplot(1,3,2); hold on;
plot(nanmean(trials_per_sess,1), 'k', 'linewidth',2);
plot(trials_per_sess', 'Color',[0.5,0.5,0.5,0.5],'linewidth',0.5)
plot(trials_per_sess(groupID==1,:)', 'Color', [.9,.2,.2,.3],'linewidth',0.5);
ylabel('Trials'); xlabel('Session');

ax3=subplot(1,3,3); hold on;
plot(nanmean(trial_time_sess,1), 'k', 'linewidth',2);
plot(trial_time_sess', 'Color',[0.5,0.5,0.5,0.5],'linewidth',0.5)
plot(trial_time_sess(groupID==1,:)', 'Color', [.9,.2,.2,.3],'linewidth',0.5);
ylabel('Trial time (s)'); xlabel('Session')

linkaxes([ax1,ax2,ax3],'x');

% plot trial times for accuracies
% figure; hold on;
% colors_use = 'krbg';
% for m = 1:length(mice)
%     subplot(2,2,m); hold on;
%     plot(acc_given_ITI_hit(m,:), 'Color', colors_use(m));
%     plot(acc_given_ITI_miss(m,:), 'Color',[0.5,0.5,0.5]);
%     legend('Hit','Miss'); 
%     xlabel('Session'); ylabel('ITI (s)')
% end

%% Temp mapke a bar plot for these metrics
nDays = 2;
acc = nanmean(accuracy_per_sess(:, end + (-nDays:0)),2);
nTrials = nanmean(trials_per_sess(:, end + (-nDays:0)),2);
trial_time = nanmean(trial_time_sess(:, end + (-nDays:0)),2);

figure;
subplot(1,3,1); hold on;
bar(mean(acc)); 
plot(ones(1,sum(groupID==0)), acc(groupID==0),'ok');
plot(ones(1,sum(groupID==1)), acc(groupID==1),'or');
ylabel('Fracion correct'); ylim([0,1]); xticks([])
subplot(1,3,2); hold on;
bar(mean(nTrials)); 
plot(ones(1,sum(groupID==0)), nTrials(groupID==0),'ok');
plot(ones(1,sum(groupID==1)), nTrials(groupID==1),'or');
ylabel('Trials / session');xticks([])
subplot(1,3,3); hold on;
bar(mean(trial_time)); 
plot(ones(1,sum(groupID==0)), trial_time(groupID==0),'ok');
plot(ones(1,sum(groupID==1)), trial_time(groupID==1),'or');
ylabel('Time (s) / trial');xticks([])


%% Temp simple learning metrics
%** rework this one, but it should print out when the animal learned?
learn_sess = []; learn_trials = [];
asymptote_acc = [];

for m = 1:nMice
    acc = accuracy_cell{m};
    badid = isnan(acc);
    acc(badid) = [];
    %acc = medfilt1(acc,4,'omitnan');
    acc = movmean(acc, [3,0], 'omitnan');
    nTrials = trials_cell{m};
    nTrials(badid) = [];
    id = find(acc > 0.6);
    % filter first trials out 
    if isempty(id)
        id = 0;
    end
    if id(1)==1
        consecidx = diff(id);
        idconsec = find(consecidx>1);
        id(1:idconsec(1)) = [];
    end
    try
        learn_sess(m) = id(1)+1; % +1 filts for medfilt
        learn_trials(m) = sum(nTrials(1:id(1)+1));
        asymptote_acc(m) = acc(end); % this has been filtered so ok
    catch
        learn_sess(m) = inf; % +1 filts for medfilt
        learn_trials(m) = inf;
        asymptote_acc(m) = inf; % this has been filtered so ok
    end
end

% [mean(learn_sess), std(learn_sess)]
% [mean(learn_trials), std(learn_trials)]
% [mean(asymptote_acc), std(asymptote_acc)]

%% pull social learning



accuracy_per_sess = {};
trials_per_sess = {}; % or trials per min? I think this is better...
trial_time_sess = {}; % this might be using stem time, not door open time...
sesstime_all = {};
mouse_name_save = {};
protocol_use = {};


% save
nMice = length(miceSocial)*2;
ITI_all = cell(nMice,1); acc_prev = cell(nMice,1);

use_asymmetry = 0;

for m = 1:length(miceSocial)
    %Prob_use = prob_use_mouse(m);
    mouse_name = miceSocial{m};
    datapath = fullfile(data_path, mouse_name, 'mouseBEHstruct.mat');
    load(datapath);
    
    for s = 1:length(info) % (length(info)-sess_use) : length(info)%
        % skip if social minus?
        if ~use_social_minus && info(s).protocolNum==3
            continue;
        end

        % load data
        LTProb = info(s).LTProb;
        choice = info(s).choice;
        reward = info(s).reward;
        badid = cellfun(@isempty, choice) | cellfun(@isempty, reward) | ...
            cellfun(@length, choice) ~= cellfun(@length, reward);
        LTProb(badid) = []; choice(badid) = []; reward(badid) = [];
        % if info(s).protocolNum == 3
        %     LT2Prob = info(s).LT2Prob;
        %     LT2Prob(badid) = [];
        % end
        if isnan(info(s).LT2Prob(1)) && length(info(s).LT2Prob)==1
            LT2Prob = LTProb;
        else
            LT2Prob = info(s).LT2Prob;
            LT2Prob(badid) = [];
        end
        
        % skip if asymmetrical probs
        if ~use_asymmetry && ~isequal(LT2Prob, LTProb) && ~use_social_minus; continue; end

        % temp for box with out color sensor
        if isfield(info(s), 'mouseID')
            mouseID = info(s).mouseID; mouseID(badid) = [];
        else
            disp('using rig 1!');
            mouseID = cellfun(@(v) (1:length(v)) - 1, choice, 'un', 0);
            % new version where just assume mice stick to top or bottom
            mouseID = choice;
            for trial = 1:length(choice)
                for i = 1:length(choice{trial})
                    if ismember(choice{trial}(i), [0,2])
                        mouseID{trial}(i) = 0;
                    elseif ismember(choice{trial}(i), [1,3])
                        mouseID{trial}(i) = 1;
                    else
                        disp('should never be here')
                    end
                end
            end
        end
            
        start_time = info(s).start_time; start_time(badid) = [];

        % skip if something bad?
        if sum(badid)/length(badid) > 0.5
            continue;
        end
        
        % format data
        choice1 = NaN(1,length(choice));  % left choice is 1, right is -1!
        reward1 = NaN(1,length(choice)); 
        choice2 = NaN(1,length(choice)); 
        reward2 = NaN(1,length(choice)); 
        outcomevals = [-1, 1];
        for trial = 1:length(choice)
            for i = 1:length(choice{trial})
                 if mouseID{trial}(i) == 0
                    choice1(trial) = choice{trial}(i);
                    reward1(trial) = reward{trial}(i);
                elseif mouseID{trial}(i) == 1
                    choice2(trial) = choice{trial}(i);
                    reward2(trial) = reward{trial}(i);
                end
            end
        end


        % save stuff

        accuracy1 = (choice1>=2 & LTProb > 50) | (choice1<=1 & LTProb<50);
        accuracy2 = (choice2>=2 & LTProb > 50) | (choice2<=1 & LTProb<50);
        if info(s).protocolNum  == 3
            accuracy2 = (choice2>=2 & LT2Prob > 50) | ...
                (choice2<=1 & LT2Prob < 50);
        end

        mm = (2*m)-1;

        % accuracy over sessions
        accuracy_per_sess{mm}(s) = mean(accuracy1);
        accuracy_per_sess{mm+1}(s) = mean(accuracy2);
        
        % number of trials over sessions
        trials_per_sess{mm}(s) = length(accuracy1);
        trials_per_sess{mm+1}(s) = length(accuracy2);
    
        % % trial time over sessions
        % trial_time = choice_time - start_time;
        % trial_time_sess{m}(s) = mean(trial_time) / 1000;

        sesstime_all{mm}(s) = info(s).sessionTime;

        protocol_use{mm}(s) = info(s).protocolNum;
        protocol_use{mm+1}(s) = info(s).protocolNum;
    
        % uh add in to get id of trianing trials?
        % if ~any(LTProb==Prob_use) || length(unique(LTProb))==1
        %     accuracy_per_sess{m}(s) = nan;
        %     trials_per_sess{m}(s) = nan;
        %     trial_time_sess{m}(s) = nan;
        % end

        % save ITI_all and accuracy prev?
        ITI = (start_time(2:end) - start_time(1:end-1))/1000;
        accuracy1_prev = accuracy1(2:end);
        accuracy2_prev = accuracy2(2:end);

        ITI_all{mm} = [ITI_all{m} , ITI];
        ITI_all{mm+1} = [ITI_all{mm+1}, ITI];
        acc_prev{mm} = [acc_prev{mm}, accuracy1_prev];
        acc_prev{mm+1} = [acc_prev{mm+1}, accuracy2_prev];
        
        mouse_name_split = strsplit(mouse_name, {username,'-'});
        mouse_name_save{mm} = [username mouse_name_split{2}];
        mouse_name_save{mm+1} = [username mouse_name_split{3}];
    
    end
end
% rearrange
for m = 1:(2*length(miceSocial))
    protocol_use{m}(isnan(accuracy_per_sess{m})) = [];
    accuracy_per_sess{m}(isnan(accuracy_per_sess{m})) = [];
    trials_per_sess{m}(isnan(trials_per_sess{m})) = [];
   
end

if use_social_minus
% pad specially for social minus
% so this pads at the transition point until the largest number of trials?
nn = min(cellfun(@(v) sum(v==3), protocol_use)); % how many sessions to take?
startIdx = cellfun(@(v) v(1), cellfun(@(v) find(v==3), protocol_use, 'un', false));
% Function to extract the values
extractFunction = @(vec, startIdx) vec((startIdx-nn+1): (startIdx+nn-1));
% Use cellfun with an additional argument for the start indices
accuracy_per_sess = cellfun(extractFunction, accuracy_per_sess, num2cell(startIdx), 'UniformOutput', false);
accuracy_per_sess = cell2mat(accuracy_per_sess');

trials_per_sess = cellfun(extractFunction, trials_per_sess, num2cell(startIdx), 'UniformOutput', false);
trials_per_sess = cell2mat(trials_per_sess');
else
% do pad?
% this pads from the beginning to the end, so looks over learning
nn = min(cellfun(@length, accuracy_per_sess));
accuracy_per_sess = cellfun(@(v) v(1:nn), accuracy_per_sess,'un',0);
accuracy_per_sess = cell2mat(accuracy_per_sess');
nn = min(cellfun(@length, trials_per_sess));
trials_per_sess = cellfun(@(v) v(1:nn), trials_per_sess,'un',0);
trials_per_sess = cell2mat(trials_per_sess');
protocol_use = cellfun(@(v) v(1:nn), protocol_use, 'un', 0);
protocol_use = cell2mat(protocol_use');

end

accuracy_per_sess_social = accuracy_per_sess;
trials_per_sess_social = trials_per_sess / 90; % normalize by session length?



%% plot both accuracies

figure; hold on;
colors_use = 'krbgymkrbgymkrbgym';
% this should crop non social sessions at start of social sessions
nback = min(size(accuracy_per_sess_social,2), size(accuracy_per_sess_nonsocial,2)) - 1;


nMice = size(accuracy_per_sess_social,1);
%nback = 6;
for m = 1:size(accuracy_per_sess_social,1)
   subplot(2,nMice/2,m); hold on;
   ax1=plot(accuracy_per_sess_nonsocial(m, end-nback:end), 'Color', [.7,.7,.7]);
   ax2=plot(accuracy_per_sess_social(m, end-nback:end), 'Color', colors_use(m));
   xlabel('Session'); ylabel('Accuracy'); %ylim([0.1,1])
   %title(mice{mod(m-1,nMice/2)+1});
   title(mouse_name_save{m});
   ylim([0,1]);
end
legend([ax1,ax2],{'non-social','social'})

% Plot all accuracies on one plot
figure; hold  on;
acc_val = [];
for m = 1:size(accuracy_per_sess_nonsocial,1)
    plot(accuracy_per_sess_nonsocial(m, end-nback:end), 'Color', [.7,.7,.7],'linewidth',0.5);
    acc_val(m,:) = accuracy_per_sess_nonsocial(m, end-nback:end);
end
ax1=plot(nanmean(acc_val,1), 'k','linewidth',2.5);
acc_val = [];
for m = 1:size(accuracy_per_sess_social,1)
    plot(accuracy_per_sess_social(m, end-nback:end), 'Color', [1,.3,.3],'linewidth',0.5);
    acc_val(m,:) = accuracy_per_sess_social(m, end-nback:end);
end
ax2=plot(nanmean(accuracy_per_sess_social,1),'r','linewidth',2.5);
xlabel('Session'); ylabel('Accuracy');

if use_social_minus
    plot(nn*[1,1], ylim, 'k--');
end

legend([ax1,ax2],{'Non social','Social'})