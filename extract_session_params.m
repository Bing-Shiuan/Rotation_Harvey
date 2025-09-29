function [choice, reward, choice_time, choice2, reward2, choice_time_2, LTProb, LTProb2] ...
    = extract_session_params(info, sess_use, RLflag)
if nargin < 2
    sess_use = 1;
end
if nargin < 3
    RLflag = 0;
end

info = info(sess_use);

if length(info)==1
    if contains(info.mouseName, '-')
        % load data
        choice = info.choice; choice_time = info.choice_time;
        reward = info.reward;
        choice_time_inferred = info.choice_time_inferred;
        badid = cellfun(@isempty, choice) | cellfun(@isempty, reward) | ...
            cellfun(@length, choice) ~= cellfun(@length, reward);
        choice(badid) = []; reward(badid) = [];
        mouseID = info.mouseID; mouseID(badid) = [];
        choice_time(badid) = []; choice_time_inferred(badid) = [];
        
        LTProb = info.LTProb; LTProb(badid) = [];
        LTProb2 = info.LT2Prob; LTProb2(badid) = [];
        
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
        if RLflag
            outcomevals = [1, 2];
        end
        for trial = 1:length(choice)
            for i = 1:length(choice{trial})
                if mouseID{trial}(i) == 0
                    % set so assings 1 or -1, and unsassigned is auto 0?
                    choice1(trial) = outcomevals( (double(choice{trial}(i)) > 1)+1 );
                    reward1(trial) = outcomevals( reward{trial}(i)+1 );
                    choice_time_1(trial) = choice_time{trial}(i);
                    choice_time_1_inferred(trial) = choice_time_inferred{trial}(i);
                elseif mouseID{trial}(i) == 1
                    choice2(trial) = outcomevals( (double(choice{trial}(i)) > 1)+1 );
                    reward2(trial) = outcomevals( reward{trial}(i)+1 );
                    choice_time_2(trial) = choice_time{trial}(i);
                    choice_time_2_inferred(trial) = choice_time_inferred{trial}(i);
                end
            end
        end
        % fix!
        if RLflag
            reward1 = reward1-1;
            reward2 = reward2-1;
        end
        % filter trials with one choice
        badid = choice1==0 | choice2==0;
        choice1(badid) = []; choice2(badid) = [];
        reward1(badid) = []; reward2(badid) = [];
        choice_time_1(badid) = []; choice_time_2(badid) = [];
        choice_time_1_inferred(badid) = []; choice_time_2_inferred(badid) = [];
        LTProb(badid) = []; % whoops this was missing
    
        % assign
        choice_time = choice_time_1_inferred;
        choice = choice1; reward = reward1;
        choice_time_2 = choice_time_2_inferred;
        
    else
        % I thought I fixed this to account for double choices...
        %choice_time = [info.choice_time{:}];
        %choice = double([info.choice{:}]>=2); 
        %reward = [info.reward{:}]; 
        
        %
        choice_time = info.choice_time;
        choice = info.choice;
        reward = info.reward;
        
        % find bad trials?
        badid = cellfun(@isempty, choice) | cellfun(@isempty, reward);
        choice_time(badid) = [];
        choice(badid) = [];
        reward(badid) = [];
        choice = cellfun(@(v) v(1), choice); choice = double(choice>=2);
        reward = cellfun(@(v) v(1), reward);
        choice_time = cellfun(@(v) v(1), choice_time);

        LTProb = info.LTProb; 
        LTProb(badid) = [];
        LTProb2 = LTProb;
        
        if RLflag
            choice = choice+1;
        else
            reward(reward==0)=-1;
            choice(choice==0)=-1;
        end
        
        % want to account for choice time inferred!
        choice_time = choice_time - 254; % it takes 254ms to play the buzzer, so they lick earlier
    
        % temp vars
        choice2 = nan(size(choice));
        reward2 = nan(size(reward));
        choice_time_2 = nan(size(choice_time));
    
        % well lets assign them to the same as choice 1 just to make things not
        % crash
        choice2 = choice;
        reward2 = reward;
        choice_time_2 = choice_time;
    end
elseif length(info) > 1
    % extract same data in a loop
    disp('todo!')
end



end