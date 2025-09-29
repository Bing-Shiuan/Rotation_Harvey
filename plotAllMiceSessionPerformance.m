%% Plot all session perforamnce
% generates plots of the current days mice performance

data_path = 'Z:\HarveyLab\Tier1\Kevin\Videos\';
addpath(genpath('../'))
addpath('Z:\HarveyLab\Tier1\Kevin\Analysis\20250718_backup_Cindys_PC\Utilities')



miceNonsocial = {'KM49','KM50'};
miceSocial = {'KM49-50'};


%% performance per mouse
initSeed = [];% 
for m = 1:length(miceNonsocial)

mouse_name = miceNonsocial{m};
days_back = 0;
info = createBEHstruct_nonsocial(mouse_name);

if days_back>0
    date_use = char(datetime('now') - days(days_back), 'yyMMdd');
    fileuse = dir(fullfile(data_path, mouse_name, [date_use '*.dat']));
    info = createBEHstruct_nonsocial(mouse_name, fileuse);
end

% Plot 
% plot of average reward rate, reward probs switch, and which side was
% chosen
LTProb = info.LTProb;
LBProb = info.LBProb;
choice = info.choice;
reward = info.reward;

% find bad trials?
badid = cellfun(@isempty, choice) | cellfun(@isempty, reward);
LTProb(badid) = [];
LBProb(badid) = [];
choice(badid) = [];
reward(badid) = [];
choice = cellfun(@(v) v(1), choice);
reward = cellfun(@(v) v(1), reward);

figure; hold on;
title([info.mouseName ' - ' datestr(info.sessionTime)])
plot(LTProb/100, 'k'); % chose top randomly

% Plot hits
RTchoice = find(choice==0 & reward);
RBchoice = find(choice==1 & reward);
LTchoice = find(choice==2 & reward);
LBchoice = find(choice==3 & reward);
plot(LTchoice, 1.07*ones(1,length(LTchoice)), '^', 'MarkerSize',3,'MarkerFaceColor',[0.9290 0.6940 0.1250],'Color',[0.9290 0.6940 0.1250]);
plot(LBchoice, 1.04*ones(1,length(LBchoice)), 'v', 'MarkerSize',3,'MarkerFaceColor',[0.9290 0.6940 0.1250],'Color',[0.9290 0.6940 0.1250]);
plot(RTchoice, -.04*ones(1,length(RTchoice)), '^', 'MarkerSize',3,'MarkerFaceColor',[0.9290 0.6940 0.1250],'Color',[0.9290 0.6940 0.1250]);
plot(RBchoice, -.07*ones(1,length(RBchoice)), 'v', 'MarkerSize',3,'MarkerFaceColor',[0.9290 0.6940 0.1250],'Color',[0.9290 0.6940 0.1250]);

% Plot misses
RTchoice = find(choice==0 & ~reward);
RBchoice = find(choice==1 & ~reward);
LTchoice = find(choice==2 & ~reward);
LBchoice = find(choice==3 & ~reward);
plot(LTchoice, 1.07*ones(1,length(LTchoice)), '^', 'MarkerSize',3,'Color',[.4,.4,.4]);
plot(LBchoice, 1.04*ones(1,length(LBchoice)), 'v', 'MarkerSize',3,'Color',[.4,.4,.4]);
plot(RTchoice, -.04*ones(1,length(RTchoice)), '^', 'MarkerSize',3,'Color',[.4,.4,.4]);
plot(RBchoice, -.07*ones(1,length(RBchoice)), 'v', 'MarkerSize',3,'Color',[.4,.4,.4]);


ylim([-.2,1.2]);
ylabel('Prob left'); xlabel('Trial');

% plot accuracy
%accuracy = (choice>=2 & LTProb > 50) | (choice<=1 & LTProb<50);
%acc = movmean(accuracy, [20, 0]);
%hold on; plot(acc,'r');

% plot smoothed prob choose left
pleft = choice>=2;
pleft = movmean(pleft, [5,5]);
hold on;plot(pleft,'r');

try
    initSeed(m) = info.randomSeedInit;
catch
    initSeed(m) = -1;
end

end

%% Plot social mice also?


    
for m = 1:length(miceSocial)
    
date_use = '250925';
    fileuse = dir(fullfile(data_path, miceSocial{m}, [date_use '*.dat']));
    info = createBEHstruct_social(miceSocial{m},fileuse);
% make the same choice plots 
% unfortunately, don't know which mouse is which...


doplottrialtime = 0;
LTProb = info.LTProb;
choice = info.choice;
reward = info.reward;
badid = cellfun(@isempty, choice) | cellfun(@isempty, reward) | ...
    cellfun(@length, choice) ~= cellfun(@length, reward);
LTProb(badid) = []; choice(badid) = []; reward(badid) = [];
trial_time = info.start_time;
trial_time(badid) = [];
choice_time = info.choice_time_inferred; % ??
choice_time(badid) = [];

try
    LT2Prob = info.LT2Prob;
    LT2Prob(badid) = [];
catch
    LT2Prob = NaN;
end

if isfield(info, 'mouseID')
    mouseID = info.mouseID;
    mouseID = info.mouseID; mouseID(badid) = [];
    % more badid?
    badid = find(cellfun(@isempty, mouseID));
    mouseID(badid) = [];
    trial_time(badid) = [];
    choice_time(badid) = [];
    choice(badid) = [];
    reward(badid) = [];
    LTProb(badid) = [];
% temp for box without color sensor
else
    mouseID = cellfun(@(v) (1:length(v)) - 1, choice, 'un', 0);
    % assume upper is mouse 1 and lower is mouse 2?
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

% - note, making up mouse id for now! 
% pull choice mouse 1
choice1 = []; reward1 = []; trial1 = []; 
choice2 = []; reward2 = []; trial2 = []; 
acc_all = []; acc1 = []; acc2 = [];
competition_trial = zeros(1,length(choice)); 
choice_time1 = []; choice_time2 = [];
for trial = 1:length(choice)
    acc = 0;
    for i = 1:length(choice{trial})
        if mouseID{trial}(i) == 0
            choice1(end+1) = choice{trial}(i);
            reward1(end+1) = reward{trial}(i);
            trial1(end+1) = trial;
            acc1(end+1) =  double( (choice{trial}(i)>=2 && LTProb(trial)>50) || (choice{trial}(i)<=1 && LTProb(trial)<50) );
            acc = acc + double( (choice{trial}(i)>=2 && LTProb(trial)>50) || (choice{trial}(i)<=1 && LTProb(trial)<50) );
            choice_time1(end+1) = choice_time{trial}(i);
        elseif mouseID{trial}(i) == 1
            choice2(end+1) = choice{trial}(i);
            reward2(end+1) = reward{trial}(i);
            trial2(end+1) = trial;
            if ~all(isnan(LT2Prob))
                acc2(end+1) = double( (choice{trial}(i)>=2 && LT2Prob(trial)>50) || (choice{trial}(i)<=1 && LT2Prob(trial)<50) );
                acc = acc + double( (choice{trial}(i)>=2 && LT2Prob(trial)>50) || (choice{trial}(i)<=1 && LT2Prob(trial)<50) );
            else
            acc2(end+1) = double( (choice{trial}(i)>=2 && LTProb(trial)>50) || (choice{trial}(i)<=1 && LTProb(trial)<50) );
            acc = acc + double( (choice{trial}(i)>=2 && LTProb(trial)>50) || (choice{trial}(i)<=1 && LTProb(trial)<50) );
            end
            choice_time2(end+1) = choice_time{trial}(i);
        end
    end
    if length(choice{trial})>1 && length(unique(choice{trial}))==1
        competition_trial(trial) = unique(choice{trial}) + 1;
    end
    acc = acc / length(choice{trial});
    acc_all(trial) = acc;
end

% plot

figure; hold on;
title([info.mouseName ' - ' datestr(info.sessionTime)])
if ~all(isnan(LT2Prob)) && ~isequal(LT2Prob, LTProb)
    plot(LTProb/100, 'Color',[0.9290 0.6940 0.1250])
    plot(LT2Prob/100, 'Color',[0.2431 0.5882 0.3176] )
else
    plot(LTProb/100, 'k'); % chose top randomly
end

% plot mouse 1
RhitID = choice1<2 & reward1;
LhitID = choice1>1 & reward1;
RmissID = choice1<2 & ~reward1;
LmissID = choice1>1 & ~reward1;
plot(trial1(LhitID), 1.04*ones(1,sum(LhitID)), '<', 'MarkerSize',3,'MarkerFaceColor',[0.9290 0.6940 0.1250],'Color',[0.9290 0.6940 0.1250]);
plot(trial1(RhitID), -.04*ones(1,sum(RhitID)), '>', 'MarkerSize',3,'MarkerFaceColor',[0.9290 0.6940 0.1250],'Color',[0.9290 0.6940 0.1250]);
plot(trial1(LmissID), 1.04*ones(1,sum(LmissID)), '<', 'MarkerSize',3,'Color',[.4,.4,.4]);
plot(trial1(RmissID), -.04*ones(1,sum(RmissID)), '>', 'MarkerSize',3,'Color',[.4,.4,.4]);

% plot mouse 2
RhitID = choice2<2 & reward2;
LhitID = choice2>1 & reward2;
RmissID = choice2<2 & ~reward2;
LmissID = choice2>1 & ~reward2;
plot(trial2(LhitID), 1.07*ones(1,sum(LhitID)), '<', 'MarkerSize',3,'MarkerFaceColor',[0.2431 0.5882 0.3176],'Color',[0.2431 0.5882 0.3176]);
plot(trial2(RhitID), -.07*ones(1,sum(RhitID)), '>', 'MarkerSize',3,'MarkerFaceColor',[0.2431 0.5882 0.3176],'Color',[0.2431 0.5882 0.3176]);
plot(trial2(LmissID), 1.07*ones(1,sum(LmissID)), '<', 'MarkerSize',3,'Color',[.4,.4,.4]);
plot(trial2(RmissID), -.07*ones(1,sum(RmissID)), '>', 'MarkerSize',3,'Color',[.4,.4,.4]);

% plot group accuracy?
hold on;
% plot(movmean(acc_all, [20,0]), 'r');
ax1=plot(trial1, movmean(acc1, [20,0]), 'Color',[0.9290 0.6940 0.1250]);
ax2=plot(trial2, movmean(acc2, [20,0]), 'Color',[0.2431 0.5882 0.3176]);

% plot competition?
id1 = find(competition_trial>0 & competition_trial<3);
plot(id1, -0.02*ones(1,length(id1)), 'k*')
id2 = find(competition_trial>2);
plot(id2,  1.11*ones(1,length(id2)), 'k*')

% % plot trial time
if doplottrialtime
yyaxis right;
trial_time = [0, trial_time];
plot(diff(trial_time)/1000)
ylabel('Trial time (s)')
end

legend([ax1,ax2],{'Mouse 1','Mouse 2'})


end



%% additionally plot swap rate all mice non social
plotSwapRate = 1;

if plotSwapRate

acc_mouse = {};
ITI_mouse = {};
prob_mouse = {};

prob_use_mice = [70,70,80,80];
mice = {'KM16','KM17','KM18','KM19'};
mice = {'KM16','KM19'};prob_use_mice = [70,80];

mice = {'KM10','KM11','KM12','KM13','KM16','KM19'};
prob_use_mice = [80,80,80,80,70,80];

%mice = {'KM16','KM19'}; prob_use_mice = [70,80];
%mice = {'KM20','KM21','KM22','KM23'};
%prob_use_mice = [80,80,80,80];

mice = {'KM24','KM25','KM26','KM27','KM28','KM29','KM30','KM31'};
prob_use_mice = [80,80,80,80,80,80,80,80];

mice = {'KM20','KM21','KM22','KM23','KM24','KM25','KM26','KM27','KM28','KM29','KM30','KM31'};
% mice = mouse_name_save; 
% mice = miceNonsocial
prob_use_mice = 80*ones(size(mice));

for m = 1:length(mice)

mouse_name = mice{m};
data_path = 'Z:\HarveyLab\Tier1\Kevin\Videos\';
%data_path = 'C:\Users\kgcmi\OneDrive\Documents\Postdoc\Social-T-maze\Analysis\data';

Nswaps = 30; % average over prev 20 swaps
tback = 10; % trials to look back
tforward = 20; % triasl to look forward
prob_use = 80; % probability to look at 
prob_use = prob_use_mice(m); % for the new cohort

% load all mouse beh data
mouse_file = fullfile(data_path, mouse_name, 'mouseBEHstruct.mat');
if exist(mouse_file)
    load(mouse_file)
else
    mouse_files = dir(fullfile(data_path, mouse_name, '*.dat'));
    mouse_files([mouse_files.bytes]<100000) = [];
    
    for sess = 1:length(mouse_files)
        mouse_path = fullfile(mouse_files(sess).folder, mouse_files(sess).name);
        info(sess) = createBEHstruct_nonsocial(mouse_name, mouse_path);
    end
end



% look at swaps
numswaps = [];
acc = []; % NaN(tforward+tback+1, Nswaps); % trials x swaps
ITI = []; % this is actually total trial time. 
prob_save = [];
for sess = 1:length(info)
    LTProb = info(sess).LTProb;
    choice = info(sess).choice;
    start_time = info(sess).start_time;
    reward = info(sess).reward;
    badid = cellfun(@isempty, choice) | cellfun(@isempty, reward);
    LTProb(badid) = []; choice(badid) = []; reward(badid) = [];
    choice = cellfun(@(v) v(1), choice);
    start_time(badid) = [];
    start_time(1) = []; start_time(end+1) = start_time(end); % um...
    trial_time = diff(start_time);
    trial_time(end+1) = trial_time(end); % buffer...

    % update if probs ever not how they are
    accuracy = (choice>=2 & LTProb > 50) | (choice<=1 & LTProb<50);
    
    swapid = find(abs(diff(LTProb)) == (2*prob_use - 100));
    swapid = find(abs(diff(LTProb)) > 20); % better way to find swaps
    swapid(swapid+tforward>length(accuracy))=[];

    %disp(length(swapid));

    for t = 1:length(swapid)
        acc(:, end+1) = accuracy(swapid(t) + (-tback:tforward));
        ITI(:, end+1) = trial_time(swapid(t) + (-tback:tforward));
        prob_save(:,end+1) = max(LTProb);
    end

end

acc_mouse{m} = acc;
ITI_mouse{m} = ITI;
prob_mouse{m} = prob_save;

end

acc_nonsocial = acc_mouse;
% filter by prob use mouse?
disp('Filtering to use only one probability condition!')
for m = 1:length(mice)
    badid = prob_mouse{m} ~= prob_use_mice(m);
    acc_nonsocial{m}(:,badid) = [];
end


Nswaps = 30;
figure; hold on;
acc_val = [];
for m = 1:length(mice)
    acc = acc_mouse{m};
    if prob_use_mice(m)<80; coloruse = 'r'; else coloruse = 'k'; end
    plot(-tback:tforward, mean(acc(:, end-Nswaps:end),2),coloruse,'linewidth',0.5); 
    acc_val(m,:) = mean(acc(:, end-Nswaps:end),2); 
end

plot(-tback:tforward, mean(acc_val),'k', 'linewidth',2.5)

ylabel('Accuracy'); xlabel('Trials to swap');
ylim([0,1])

% temp to look at ITI?
% Nswaps = 30;
% figure; hold on;
% ITI_val = [];
% for m = 1:length(mice)
%     ITI = ITI_mouse{m};
%     plot(-tback:tforward, mean(ITI(:, end-Nswaps:end),2),'k','linewidth',0.5); 
%     ITI_val(m,:) = mean(ITI(:, end-Nswaps:end),2); 
% end
% plot(-tback:tforward, mean(ITI_val),'k', 'linewidth',2.5)

% look only at slow and fast trials?
% - ok so mice relearn more slowly when there are long ITIs.
% - is the frequency of ITIs higher in social mice?
colors_use = 'krbgymkrbgym';
figure; 
accslow_nonsocial = {};
accfast_nonsocial = {};
for m = 1:length(mice)
    subplot(2,length(mice)/2,m); hold on;
    ITI = ITI_mouse{m};
    acc = acc_mouse{m};
    
    id = any(ITI((tback+1):end, :) > 30*1000);
    accslow = acc(:,id);
    accfast = acc(:,~id);
    accslow_nonsocial{m} = accslow;
    accfast_nonsocial{m} = accfast;
    
    plot(-tback:tforward, mean(accslow(:, end-Nswaps:end),2),'Color',[0.5,0.5,0.5]); 
    plot(-tback:tforward, mean(accfast(:, end-Nswaps:end),2),'Color', colors_use(m)); 

    xlabel('Trials to swap'); ylabel('accuracy'); title(mice{m});
    legend('> 30 seconds ITI','< 30 seconds ITI')
end

acc_nonsocial = acc_mouse;
ITI_nonsocial = ITI_mouse;

% plot early vs. late

Nswaps = 30;
acc_late = []; acc_early = [];
figure; hold on;
for m = 1:length(mice)
    acc = acc_mouse{m};
    plot(-tback:tforward, mean(acc(:, 1:Nswaps),2),'Color',[0.5,0.5,0.5],'linewidth',0.5); 
    plot(-tback:tforward, mean(acc(:, end-Nswaps:end), 2), 'Color', [0.7,0,0], 'linewidth',0.5);
    acc_late(m,:) = mean(acc(:, end-Nswaps:end),2); 
    acc_early(m,:) = mean(acc(:,1:Nswaps), 2);
end

ax1=plot(-tback:tforward, mean(acc_early),'Color',[0.5,0.5,0.5], 'linewidth',2.5);
ax2=plot(-tback:tforward, mean(acc_late) ,'Color',[0.7,0,0], 'linewidth',2.5);
ylabel('Accuracy'); xlabel('Trials to swap');
legend([ax1,ax2],{'Early','Late'})
title(['Swap rate averaged over N = ' num2str(Nswaps)])
ylim([0,1])


end
