%% Check the water drank during session

% set path to .dat files!!!
fileLoc = 'Z:\HarveyLab\Tier1\Kevin\Videos'; 

miceNonsocial = {'KM49','KM50'};
miceSocial = {'KM49-50'};

% miceSocial = {'KM51-52','KM55-56','KM57-58','KM61-62','KM66-67'};
% miceNonsocial = {'KM49','KM50','KM59','KM60','KM63','KM64', 'KM65'};


% hard coded variables
ulPerPulse = 3; % 4 ul of water per each pulse <- need to recalibrate might be like 3
WaterIntakePIN    = 51;
TrialCountPIN     = 97;

%% NON SOCIAL MICE



% Loop through mice
totalWater = zeros(size(miceNonsocial));
numTrials = zeros(size(miceNonsocial));
for m = 1:length(miceNonsocial)
    % Get filepath for most recent .dat
    filename = dir(fullfile(fileLoc, miceNonsocial{m}, '*.dat'));
    [~,idx] = max([filename.datenum]);
    filename = fullfile(filename(idx).folder, filename(idx).name);

    % Instead use the most recent date
    filename = dir(fullfile(fileLoc, miceNonsocial{m}, '*.dat'));
    names = {filename.name};
    pattern = '^\d{12}\.dat$';
    datFolders = names(~cellfun(@isempty, regexp(names, pattern)));
    timestamps = datetime(datFolders, 'InputFormat', 'yyMMddHHmmss''.dat''');
    [~, idx] = max(timestamps);
    filename = fullfile(filename(idx).folder, filename(idx).name);


    sessionTime = regexp(filename, '\d{12}', 'match');
    sessionTime = datetime(sessionTime{1}, 'InputFormat', 'yyMMddHHmmss');
    %disp(sessionTime)
    if hours(datetime('now') - sessionTime) > 12
        disp('This file was over 12 hours ago!')
    end

    % Read output file
    fileID = fopen(filename, 'r');
    dataArray = textscan(fileID, '%f%f%f', 'Delimiter',',');
    fclose(fileID);
    nn = min(cellfun(@length,dataArray));
    dataArray = cellfun(@(v) v(1:nn), dataArray, 'UniformOutput',false);
    EventType =    dataArray{:, 1};
    Value =        dataArray{:, 2};
    SessTime =     dataArray{:, 3};

    % Get water
    sessionWater  = Value(EventType == WaterIntakePIN);
    if isempty(sessionWater); sessionWater = 0; end
    totalWater(m) = sessionWater(end) * ulPerPulse; % water in microliters!
    %numTrials(m) = length(sessionWater); % this isn't right!!!! 
    TrialCount    = Value(EventType == TrialCountPIN);
    numTrials(m)  = length(TrialCount);

    disp([miceNonsocial{m} ', ' num2str(totalWater(m)/1000) ' ml, ' num2str(numTrials(m)) ' trials'])
    
end


%% run for social

for m = 1:length(miceSocial)

info = createBEHstruct_social(miceSocial{m},0);

totalRewards = [0;0];
% Loop through each entry in the arrays
for i = 1:length(info.reward)
    curreward = info.reward{i};
    curmouse = info.mouseID{i};
    if isempty(curmouse); continue; end
    if isempty(curreward); continue; end
    assert(length(curmouse)==length(unique(curmouse)))
    minn = min(length(curreward), length(curmouse));
    curreward = curreward(1:minn); curmouse = curmouse(1:minn); % hacky fix
    totalRewards(curmouse+1) = totalRewards(curmouse+1) + curreward;
    % alternate way
    % for j = 1:length(curreward)
    %     totalRewards(curmouse(j)+1) = totalRewards(curmouse(j)+1) + curreward(j);
    % end
end
totalWater = totalRewards * ulPerPulse / 1000;

if hours(datetime('now') - info.sessionTime) > 12
    disp('This file was over 12 hours ago!')
end

disp([info.mouseName, ', ', num2str(totalWater(1)) ' ml, ', ...
    num2str(totalWater(2)), ' ml, ' ...
    num2str(length(info.reward)), ' trials'])

end

%% Look at proportion of trials
if 0
load('Z:\HarveyLab\Tier1\Kevin\Videos\KM37-38\mouseBEHstruct.mat')

fracN = length(info);
for s = 1:length(info)

[choice, reward, choice_time, choice2, reward2, choice_time2, LTProb] = extract_session_params(info(s));
[uniqueChoices, ~, choiceIDX] = unique([choice; choice2]', 'rows');
minN = min(accumarray(choiceIDX,1));
unique_choices = accumarray(choiceIDX,1);
totalN = sum(unique_choices);
fracN(s) = minN / totalN;

end

p = [info(:).protocolNum];

figure; plot(fracN);
yyaxis right; plot(p); ylim([1,5])

end