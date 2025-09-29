function info = createBEHstruct_nonsocial( mouseName, filename, boxtype, skipvid)

%TODO: this function returns a lot of nans. there's some issue detecting
%the trials, or theres dropped info. Need to fix this up in the future.
%Likely an issue with some low-level coding of how a trial is defined...

info = struct;
%% Generate the behavioral structures
% Inputs: location of a experiment containing .dat files
%  - and potentially .events files to sync with the camera
%  - and potentially a file with info on what pins correspond to what ports
% Outputs: a .mat file with mouseBEH structure containing fields:
%   mouseName
%   sessionTime
%   protocolNum
%   totalWater      total number of water pulses mouse drinks
%   choice          vector of L(0,1)/R(2,3) choice for each trial 
%   reward          vector of 0/1 for whether reward was given
%   P_RT            vector of port probability for giving reward
%   P_RB
%   P_LT
%   P_LB
%   choice_time     vector of time in ms when choice was made (lick)
%   start_time      vector of time in ms when the trial started (door)
%   trial_time      vector of time in ms when mouse initiated trial (IR)
%   end_time        vector of time in ms when mouse returned
%  TODO: all times in frames
%

%% 
if contains(mouseName,'KM')
    data_path = 'Z:\HarveyLab\Tier1\Kevin\Videos';
    if contains(mouseName, {'KM57','KM58','KM59','KM60'}) 
        data_path = 'Z:\HarveyLab\Tier1\Rhyanne';
    end
elseif contains(mouseName, 'RF') 
    data_path = 'Z:\HarveyLab\Tier1\Rhyanne';
end

if nargin < 2 || strcmp(filename, 'today')
    % pull most recent file
    filename = dir(fullfile(data_path, mouseName, '*.dat'));
    [~,idx] = max([filename.datenum]);
    filename = fullfile(filename(idx).folder, filename(idx).name);
end

if isnumeric(filename)
    % assume this is days off
    days_back = filename;
    date_use = char(datetime('now') - days(days_back), 'yyMMdd');
    %data_path = 'Z:\HarveyLab\Tier1\Kevin\Videos';
    filename = dir(fullfile(data_path, mouseName, [date_use '*.dat']));
    % fix for multiple files of the same day
    % - by taking the most recent
    [~,idx] = max([filename.datenum]);
    filename = fullfile(filename(idx).folder, filename(idx).name);
end

if isstruct(filename)
    filename = fullfile(filename.folder, filename.name);
end

if nargin < 3
    boxtype = 'None';
end
if nargin< 4
    skipvid = 1;
end



%% Hardcoded orders
choice_labels = {'Right top','Right bottom','Left top','Left bottom'};

%% Pull hardware pin info

% hardcoded info
% none of these are maze specific so can use across mazes?
SessionPIN        = 99;
TrialPIN          = 90;
TrialCountPIN     = 97;
SessionStartupPIN = 100;
WaterIntakePIN    = 51;
ChoicePIN         = 91;
RewardPIN         = 92;
RTProbPIN         = 93;
RBProbPIN         = 94;
LTProbPIN         = 95;
LBProbPIN         = 96;
ProtocolPIN       = 101; % now defined as 0 for port trianing, 1 for non-social, 2 for social, 3 for anti social

% Box specific pins!
% RFID - this is pin specific... ugh
RedPIN            = 251;
GreenPIN          = 252;
BluePIN           = 253;
ClearPIN          = 254;
Dist1PIN          = 350;
Dist2PIN          = 351;
Dist3PIN          = 352;
Dist4PIN          = 353;

% BOX SPECIFIC PINS
RFIDRTPIN = 15;
RFIDLTPIN = 16;


%% Read output file
delimiter = ',';

% Format string for each line of text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f';

% Open the text file.
fileID = fopen(filename,'r');

% Read columns of data according to format string.
try
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
catch
    dataArray = cell(1,8);
end

% Close the text file.
fclose(fileID);

%fix different lengths
nn = min(cellfun(@length,dataArray));
dataArray = cellfun(@(v) v(1:nn), dataArray, 'UniformOutput',false);

%* changeme
EventType =    dataArray{:, 1};
Value =        dataArray{:, 2};
SessTime =     dataArray{:, 3};

%% sync events with video frames
% ** NOTE: THESE ARE MAZE SPECIFIC! **
SessFrame = NaN(size(SessTime));

if ~skipvid
% For maze 1
% TEENSY VIDEO PINS
cam1 = 38;
cam2 = 39;
cam3 = 40;
cam4 = 41;
cam5 = 13;
cam6 = 14;

% PICAMERA pins
pin1 = 0;
pin2 = 5;
pin3 = 6;
pin4 = 13;
pin5 = 19;
pin6 = 26; % this is not used for some cams

if strcmp(boxtype, 'T-maze2')
    % TEENSY VIDEO PINS
    cam1 = 6; cam2 = 5; cam3 = 4; cam4 = 3; cam5 = 2; cam6 = 1;
end

% pull files
filematch = strsplit(filename, '\'); 
filematch = filematch{end}(1:6);
videvents = dir(fullfile('Z:\HarveyLab\Tier1\Kevin\Videos', mouseName, [filematch '*.events']));

%if length(videvents)==1
try
videvents = fullfile(videvents.folder, videvents.name);

% read output file
delimiter = ',';
formatSpec = '%f%f%f%f';
% Open the text file.
fileID = fopen(videvents,'r');
% Read columns of data according to format string.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
% Close the text file.
fclose(fileID);

%fix different lengths
nn = min(cellfun(@length,dataArray));
dataArray = cellfun(@(v) v(1:nn), dataArray, 'UniformOutput',false);

% Save data
CamEvent =    dataArray{:, 1};
CamValue =    dataArray{:, 2};
CamFrame =    dataArray{:, 3};
CamTime  =    dataArray{:, 4};

% sync -- issue is some of the teensy events get dropped!
TeensyIDX = ismember(EventType, [cam1, cam2, cam3, cam4, cam5, cam6]);
TeensyEvent = EventType(TeensyIDX);
TeensyValue = Value(TeensyIDX);
TeensyTime  = SessTime(TeensyIDX);

% Function to sync Cam and Teensy Time
% - matches differences in timing events
% - generates a best fit line 
% - returns, SessFrame, converting SessTime to frames based on this fit
[SessFrame,bTeensyToCam, bCamToFrame, bTeensyToFrame] = ...
    syncPiTeensy(TeensyTime(TeensyValue==1), SessTime, CamTime(CamValue==1), CamFrame(CamValue==1));
%SessFrame = syncPiTeensy(TeensyTime, SessTime, CamTime, CamFrame);
catch
disp('Issue syncing frames. Check if theres an .events file...')
SessFrame = NaN(size(SessTime));
end
end

%% Save session info
protocolNum = Value(EventType == ProtocolPIN);
sessionWater  = Value(EventType == WaterIntakePIN);
totalWater = max(sessionWater) * 4; % water in microliters!
% parse date from filename
sessionTime = regexp(filename, '\d{12}', 'match');
sessionTime = datetime(sessionTime{1}, 'InputFormat', 'yyMMddHHmmss');
   

%% Water conditioning
if protocolNum == 0
    % run a script since there's no trial structure!
    % - just count choice licks, and save times, and save port IDs
    % vector of 0/1/2/3, vector of times, 
    
    % pull all water delivery times!
    idxWaterIntake = find(EventType == WaterIntakePIN);
    choice = []; reward = []; choice_time = [];
    for j = 1:length(idxWaterIntake)
        % find which port was licked right before
        idx = find(EventType(1:idxWaterIntake(j)) == ChoicePIN); idx = idx(end);
        choice(j) = Value(idx);
        reward(j) = 1;
        choice_time(j) = SessTime(idx);
    end
    
    NumTrials = length(idxWaterIntake);

    % other variables
    RTProb = ones(1,NumTrials);
    RBProb = ones(1,NumTrials);
    LTProb = ones(1,NumTrials);
    LBProb = ones(1,NumTrials);
    trial_time = [];
    end_time   = [];
    start_time = [];

end

%% Save trial info

if protocolNum == 1 || protocolNum == 2 

% pull trial start times and trial end times
TrialStartIDX = EventType == TrialPIN & Value == 1; % when mouse breaks beam
TrialStopIDX  = EventType == TrialPIN & Value == 0;
TrialBeginIDX = EventType == TrialCountPIN; %*** this might be dropping outputs!!! 
NumTrials     = Value(EventType == TrialCountPIN); % this is definitely dropping outputs...
% check pulled the right number of trials?
%assert(sum(TrialStartIDX)==length(NumTrials))



NumTrials = length(NumTrials);
TrialStart = find(TrialStartIDX);
TrialStop = find(TrialStopIDX);
TrialBegin = find(TrialBeginIDX);

% temp for social?
TrialStart = TrialBegin;
TrialStop = [TrialBegin(2:end) - 1; length(TrialBeginIDX)];

% Check! Trial Start and Trial STop may be misaligned
if length(TrialStart)+1 == length(TrialStop)
    if TrialStop(1) < TrialStart(1)
        disp('Missing a start maybe. Adding in a fake first trial that will be cut')
        TrialStart = [1; TrialStart]; % fake trial...
    end
end
if length(TrialStop)+1 == length(TrialStart)
    if all((TrialStop - TrialStart(1:end-1) )>0)
        disp('Missing end maybe. Removing final  trial')
        TrialStart(end) = [];
    end
end
%TODO fix for mouse 9 where .dat file was split into two...
if any(TrialStart > TrialStop)
    disp('Trial counter is misaligned!')
    disp('Trying to fix')
    tmp = Value(EventType==TrialPIN);
    if any(diff(tmp)==0); disp('Missed a trial signal!'); end
    if all(TrialStop < TrialStart);
        disp('shifted, adding in an extra start / stop')
        disp('This will likely result in an empty trial')
        TrialStart = [1; TrialStart];
        TrialStop = [TrialStop; length(TrialStopIDX)];
    end
end

% Simpler way for non-social, if there's only a signle choice per trial
% - this wont work due to bugs! might as well save them all.
%choice = Value(EventType == ChoicePIN);
%reward = Value(EventType == RewardPIN);

% Ver 2: the loop

% vars to save
choice = cell(1,NumTrials);
reward = cell(1,NumTrials);
RTProb = NaN(1,NumTrials);
RBProb = NaN(1,NumTrials);
LTProb = NaN(1,NumTrials);
LBProb = NaN(1,NumTrials);
choice_time = cell(1,NumTrials);
trial_time = NaN(1,NumTrials);
end_time = NaN(1,NumTrials);
start_time = NaN(1,NumTrials);
% add in frames
choice_frame = cell(1,NumTrials);
trial_frame = NaN(1,NumTrials);
end_frame = NaN(1,NumTrials);
start_frame = NaN(1,NumTrials);

if strcmp(boxtype, 'Color') || strcmp(boxtype, 'T-maze2')
    colorscan = NaN(4, NumTrials); 
    mousedist = NaN(4, NumTrials);
end
if strcmp(boxtype, 'Rfid') || strcmp(boxtype, 'T-maze1')
    LTid = NaN(1, NumTrials);
    RTid = NaN(1, NumTrials);
    mouseID = cell(1, NumTrials); % will try to infer!
    % - if two choices, one id, infer
    % - if one chocie, ignore?
end

% Loop over all trials
for trial = 1:NumTrials
    trialIDX = SessTime >= SessTime(TrialStart(trial)) & SessTime < SessTime(TrialStop(trial));

    try
    % pull data for this trial
    choice{trial}   = Value(trialIDX & EventType == ChoicePIN);
    reward{trial}   = Value(trialIDX & EventType == RewardPIN);
    RTProb(trial)   = Value(trialIDX & EventType == RTProbPIN);
    RBProb(trial)   = Value(trialIDX & EventType == RBProbPIN);
    LTProb(trial)   = Value(trialIDX & EventType == LTProbPIN);
    LBProb(trial)   = Value(trialIDX & EventType == LBProbPIN);
    choice_time{trial} = SessTime(trialIDX & EventType == ChoicePIN); % these might be wrong...
    trial_time(trial)  = SessTime(TrialStart(trial));
    end_time(trial)    = SessTime(TrialStop(trial));
    start_time(trial)  = SessTime(TrialBegin(trial));

    choice_frame{trial} = SessFrame(trialIDX & EventType == ChoicePIN);
    trial_frame(trial) = SessFrame(TrialStart(trial));
    end_frame(trial) = SessFrame(TrialStop(trial));
    start_frame(trial) = SessFrame(TrialBegin(trial));
    

    % box specific info
    % box 1 - RFID info (this is not working yet!!!)
    if strcmp(boxtype, 'RFID') || strcmp(boxtype, 'T-maze1')
        rfidRT = Value(trialIDX & EventType == RFIDRTPIN);  % this gets read on scan, and on return
        rfidLT = Value(trialIDX & EventType == RFIDLTPIN);
        rfidRT(rfidRT<2)=[]; rfidLT(rfidLT<2) = [];
        if isempty(rfidRT); rfidRT = NaN; end
        if isempty(rfidLT); rfidLT = NaN; end

        LTid(trial) = rfidLT; 
        RTid(trial) = rfidRT;
        
        % make inferences for two mice?
        idLT = find(choice{trial}>1);
        idRT = find(choice{trial}<2);
        miceID{trial} = [0, 0];
        miceID{trial}(idLT) = rfidLT;
        miceID{trial}(idRT) = rfidRT;

    % box 2 - color info
    elseif strcmp(boxtype, 'Color') || strcmp(boxtype, 'T-maze2')
        r = Value(trialIDX & EventType == RedPIN);
        g = Value(trialIDX & EventType == GreenPIN);
        b = Value(trialIDX & EventType == BluePIN);
        c = Value(trialIDX & EventType == ClearPIN);
        r = r(1); g = g(1); b = b(1); c = c(1);
        colorscan(:,trial) = [r,g,b,c];
        % also save mouse distance?
        m1 = Value(trialIDX & EventType == Dist1PIN);
        m2 = Value(trialIDX & EventType == Dist2PIN);
        m3 = Value(trialIDX & EventType == Dist3PIN);
        m4 = Value(trialIDX & EventType == Dist4PIN);
        if isempty(m1); m1 = NaN; end
        if isempty(m2); m2 = NaN; end
        if isempty(m3); m3 = NaN; end
        if isempty(m4); m4 = NaN; end
        m1 = m1(1); m2 = m2(1); m3 = m3(1); m4 = m4(1);
        mousedist(:,trial) = [m1,m2,m3,m4];
    end


    catch
        disp('some problem to fix for later')
    end
end

end

%% fix up mice ID
if strcmp(boxtype, 'RFID') || strcmp(boxtype, 'T-maze1')
    RFIDtags = unique([miceID{:}]);
    RFIDtags(isnan(RFIDtags)) = [];
    RFIDtags(RFIDtags<10) = [];
    for trial = 1:length(miceID)
        miceID{trial}(miceID{trial}==RFIDtags(1)) = 1;
        miceID{trial}(miceID{trial}==RFIDtags(2)) = 2;
    end
end

%% dynamically assign protocol num
if length(unique([LTProb, RTProb]))==1
    protocolNum=1; % no bandit swapping
else
    protocolNum=2; % bandit swapping
end

%% save into struct
info.mouseName = mouseName;
info.sessionTime = sessionTime;
info.protocolNum = protocolNum;
info.totalWater = totalWater;
info.NumTrials = NumTrials;
info.choice = choice;
info.reward = reward;
info.RTProb = RTProb;
info.RBProb = RBProb;
info.LTProb = LTProb;
info.LBProb = LBProb;
info.choice_time = choice_time;
info.trial_time = trial_time;
info.start_time = start_time;
info.choice_labels = choice_labels;
info.choice_frame = choice_frame;
info.trial_frame = trial_frame;
if strcmp(boxtype, 'Color') || strcmp(boxtype, 'T-maze2')
    info.colorscan = colorscan;
    info.mousedist = mousedist;
end

% a check!
if ~skipvid
    info.sync_check = CamFrame(CamValue == 1 & CamEvent == 0);
    info.bTeensyToCam = bTeensyToCam;
    info.bCamToFrame = bCamToFrame;
    info.bTeensyToFrame = bTeensyToFrame;
end

% for python
info.sessionTime4Python = char(sessionTime);


end



