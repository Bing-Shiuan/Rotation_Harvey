function info = createBEHstruct_social( mouseName, filename, boxtype, skipvid)

%% Function to process the social data
% - this currently only works for color, since I didn't get the RFID
% working
% output is very similar to previous mouse struct, butinclude an additional
% param that contains the ID of which mouse did wich action



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
info = struct;

if contains(mouseName,'KM')
    data_path = 'Z:\HarveyLab\Tier1\Kevin\Videos';
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
end

if isstruct(filename)
    filename = fullfile(filename(1).folder, filename(1).name);
end

if nargin < 3
     % infer box type from hard coded mouse info
    if contains(mouseName, {'10','11','12','13','16','17','18','19'})
        boxtype = 'T-maze2';
    elseif contains(mouseName,{'20','21','22','23'})
        boxtype = 'T-maze3';
    elseif contains(mouseName, {'24','25','26','27'})
        boxtype = 'T-maze4';
    elseif contains(mouseName, {'28','29','30','31'})
        boxtype = 'T-maze5';
    elseif contains(mouseName, {'32','33','34','35','36','37','38','39','40','41','42'})
        boxtype = 'T-maze-ephys';
    else
        boxtype = 'None';
    end
end

if nargin < 4
    skipvid = 1;
end

%% Hardcoded orders
choice_labels = {'Right top','Right bottom','Left top','Left bottom'};

%% Hardcoded PIN info

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
ProtocolPIN       = 101;

LT2ProbPIN = LTProbPIN+10;
RT2ProbPIN = RTProbPIN+10;


% Box specific pins!
% RFID - this is pin specific... ugh
% - note that these all won't show up all the time!
RedPIN            = 251;
GreenPIN          = 252;
BluePIN           = 253;
ClearPIN          = 254;
Dist1PIN          = 350;
Dist2PIN          = 351;
Dist3PIN          = 352;
Dist4PIN          = 353;
MouseIDPIN        = 200; % <- new
MouseRewardIDPIN  = 201; % <- newer, spit out only after reward was delivered or not
PortRewardIDPIN   = 202;

% Digital pins
if strcmp(boxtype,'T-maze2')
    LickPIN = [41,14,11,12];
elseif strcmp(boxtype, 'T-maze3')
    LickPIN = [2,5,30,27];
elseif strcmp(boxtype, 'T-maze4')
    LickPIN = [36, 39, 0, 23];
elseif strcmp(boxtype,'T-maze5')
    LickPIN = [38,40,7,5];
elseif strcmp(boxtype, 'T-maze-ephys')
    LickPIN = [36,39,21,3];
else
    LickPIN = [-1,-1,-1,-1]; % no pins
end


% RFID specific pins
RTPIN = 15; % RFID pin of RT
LTPIN = 16; % RFID pin of LT

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

%% TODO, sync with video frames
SessFrame = NaN(size(SessTime));
if ~skipvid
% For maze 2
if strcmp(boxtype, 'None') || strcmp(boxtype, 'T-maze2')
% TEENSY VIDEO PINS
cam1 = 6;
cam2 = 5;
cam3 = 4;
cam4 = 3;
cam5 = 2;
cam6 = 1;

% PICAMERA pins
pin1 = 0;
pin2 = 5;
pin3 = 6;
pin4 = 13;
pin5 = 19;
pin6 = 26; % this is not used for some cams
end

% pull video files
filematch = strsplit(filename, '\'); 
filematch = filematch{end}(1:6);
videvents = dir(fullfile('Z:\HarveyLab\Tier1\Kevin\Videos', mouseName, [filematch '*.events']));
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
[SessFrame, bTeensyToCam, bCamToFrame] = syncPiTeensy(TeensyTime(TeensyValue==1), SessTime, CamTime(CamValue==1), CamFrame(CamValue==1));
%SessFrame = syncPiTeensy(TeensyTime, SessTime, CamTime, CamFrame);

end



%% Get session info

% TOP level
protocolNum = Value(EventType == ProtocolPIN);
sessionWater  = Value(EventType == WaterIntakePIN);
totalWater = max(sessionWater) * 4; % water in microliters!
% parse date from filename
sessionTime = regexp(filename, '\d{12}', 'match');
sessionTime = datetime(sessionTime{1}, 'InputFormat', 'yyMMddHHmmss');

% fix
if isempty(protocolNum)
    protocolNum = 0;
end
if length(protocolNum)>1
    protocolNum = protocolNum(1);
end

% fix protocol num pin dynamically
% - basically figure out if in social + or social minus


%% Trial info

% pull trial start times and trial end times
TrialStartIDX = EventType == TrialPIN & Value == 1; % when mouse breaks beam
TrialStopIDX  = EventType == TrialPIN & Value == 0; % is this fixed?
TrialBeginIDX = EventType == TrialCountPIN; %*this might be dropping outputs!!! 
NumTrials     = Value(EventType == TrialCountPIN); % this is definitely dropping outputs...
% check pulled the right number of trials?
%assert(sum(TrialStartIDX)==length(NumTrials))

NumTrials = length(NumTrials);
TrialStart = find(TrialStartIDX);
TrialStop = find(TrialStopIDX);
TrialBegin = find(TrialBeginIDX);

if NumTrials > length(TrialStart)
    NumTrials = length(TrialStart);
end

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

% More aggressive fix <-too aggressive?
for trial = 1:(min([length(TrialStart), length(TrialStop)])-1)
    if trial+1 > length(TrialStart); continue; end
    dt = TrialStart(trial+1) - TrialStop(trial);
    if dt < 0; TrialStart(trial+1) = []; 
        disp(trial); trial = 1;
        disp('more aggressive fix here')
    end
    if dt > 80; TrialStop(trial) = []; trial = 1; end
end
if isempty(TrialStart) || isempty(TrialStop)
    disp('I messed up real bad now')
    figure; waitforbuttonpress;
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

if NumTrials > length(TrialStart)
    NumTrials = length(TrialStart);
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
mouseID = cell(1, NumTrials);
% time of first lick <- this is a better metric of when the choice was made
choice_time_inferred = cell(1, NumTrials);
choice_frame_inferred = cell(1, NumTrials);

%if protocolNum==3
LTProb2 = NaN(1,NumTrials); LBProb2 = NaN(1,NumTrials);
RTProb2 = NaN(1,NumTrials); RTProb2 = NaN(1,NumTrials);
%end

if strcmp(boxtype, 'Color') || strcmp(boxtype, 'T-maze2')
    colorscan = NaN(4, NumTrials); 
    mousedist = NaN(4, NumTrials);
end

% Loop over all trials
for trial = 1:NumTrials
    trialIDX = SessTime >= SessTime(TrialStart(trial)) & SessTime < SessTime(TrialStop(trial));

    try
    % want to pull the port that was licked right when mouseID is deliverd?

    % pull data for this trial
    choice{trial}   = Value(trialIDX & EventType == ChoicePIN);
    reward{trial}   = Value(trialIDX & EventType == RewardPIN);
    RTProb(trial)   = Value(trialIDX & EventType == RTProbPIN);
    RBProb(trial)   = Value(trialIDX & EventType == RBProbPIN);
    LTProb(trial)   = Value(trialIDX & EventType == LTProbPIN);
    LBProb(trial)   = Value(trialIDX & EventType == LBProbPIN);
    %if protocolNum==3
    try
        LTProb2(trial) = Value(trialIDX & EventType == LT2ProbPIN);
        RTProb2(trial) = Value(trialIDX & EventType == RT2ProbPIN);
    catch
        assert(isempty(find(EventType==LT2ProbPIN)),'something wrong with port probs pins');
        LTProb2(trial) = LTProb(trial);
        RTProb2(trial) = RTProb(trial);
    end
    choice_time{trial} = SessTime(trialIDX & EventType == ChoicePIN);
    trial_time(trial)  = SessTime(TrialStart(trial));
    end_time(trial)    = SessTime(TrialStop(trial));
    start_time(trial)  = SessTime(TrialBegin(trial));

    

    % old! all mice ids
    mouseID{trial}     = Value(trialIDX & EventType  == MouseIDPIN);

    % METHOD 1: just take the first unique mouseID and count that
    [tmp,bid] = unique(mouseID{trial}, 'stable');
    bid(tmp==-1) = [];
    mouseID{trial} = mouseID{trial}(bid);
    choice{trial} = choice{trial}(bid);
    choice_time{trial} = choice_time{trial}(bid);


    % METHOD 2: take the info just before the reward was delivered

    % METHOD 3: added in some new pins (201, and 202) which will hopefully
    % be more accurate
     % - delivered right after reward, and should be delivered only twice a
     % session!!!! 
    mouseID{trial} = Value(trialIDX & EventType  == MouseRewardIDPIN);
    choice{trial}  = Value(trialIDX & EventType == PortRewardIDPIN);
    choice_time{trial} = SessTime(trialIDX & EventType == PortRewardIDPIN);

    if any(mouseID{trial} > 1)
        disp('Incorrect mouseRewardIDpin');
        disp('trying a fix');
        badid = find(mouseID{trial}>1);
        goodid = find(mouseID{trial}==0 | mouseID{trial}==1);
        mouseID{trial}(badid) = 1-mouseID{trial}(goodid); % ugh
    end

    % add in first lick time
    % V1: assume mouse to trigger port is first mouse to lick
    for ccc = 1:length(choice{trial})
        trialLickTimes  = SessTime(trialIDX & EventType == LickPIN(choice{trial}(ccc)+1));
        % V1: take first lick (assuming mouse is first to trigger)
        if ~isempty(trialLickTimes); choice_time_inferred{trial}(ccc) = trialLickTimes(1);
        else choice_time_inferred{trial}(ccc) = NaN;
        end
        % V2: take lick nearest the reward time <-but a lot of licks can
        % happen in the meantime
    end

    % debug
    if length(choice{trial}) ~= length(reward{trial})
        disp('a reward or choice info got dropped...')
        %reward{trial} = reward{trial}(1:length(choice{trial}));
    end

    if length(mouseID{trial}) > 2
        mouseID{trial} = unique(mouseID{trial},'stable');
    end

    if length(mouseID{trial})~=length(choice{trial}) && length(mouseID{trial})==1
        disp('attempting to fix a dropped mouseID!')
        % Find when reported, and assuming missing one is further from
        missingMouseID = setdiff(0:1, mouseID{trial});
        tmpIDPort = find(trialIDX & EventType == PortRewardIDPIN);
        tmpIDMouse = find(trialIDX & EventType == MouseRewardIDPIN);
        if abs((tmpIDPort(1) - tmpIDMouse)) < abs((tmpIDPort(2) - tmpIDMouse))
            mouseID{trial} = [mouseID{trial} ; missingMouseID];
        else
            mouseID{trial} = [missingMouseID ; mouseID{trial}];
        end
    end

    if ~skipvid
        %*** todo, filter these by the reward pins?!?!?!?!?!
    choice_frame{trial} = SessFrame(trialIDX & EventType == PortRewardIDPIN); % this was choice pin
    trial_frame(trial) = SessFrame(TrialStart(trial));
    end_frame(trial) = SessFrame(TrialStop(trial));
    start_frame(trial) = SessFrame(TrialBegin(trial));
    end

    % box specific info
    % box 1 - RFID info (this is not working yet!!!)
    % box 2 - color info
    r = Value(trialIDX & EventType == RedPIN);
    g = Value(trialIDX & EventType == GreenPIN);
    b = Value(trialIDX & EventType == BluePIN);
    c = Value(trialIDX & EventType == ClearPIN);
    if isempty(r); r = NaN; end
    if isempty(g); g = NaN; end
    if isempty(b); b = NaN; end
    if isempty(c); c = NaN; end
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


    catch
        disp('some problem to fix for later')
    end
end

%% Dynamically assign protocol num...
if isequaln(LTProb2, LTProb)
    protocolNum = 2; % social plus
else
    tmpa = sign(diff(LTProb));
    tmpb = sign(diff(LTProb2));
    if isequaln(tmpa, tmpb)
        protocolNum = 2; % social plus with asymmetry
    else
        if isequal(find(tmpa~=0), find(tmpb~=0))
            protocolNum = 3; % social minus!
        else
            protocolNum = 4; % social neutral!
        end
    end
end

%% save

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
%if protocolNum == 3
info.LT2Prob = LTProb2;
info.RT2Prob = RTProb2;
%else
%    info.LT2Prob = NaN;
%    info.RT2Prob = NaN;
%end
info.choice_time = choice_time;
info.choice_time_inferred = choice_time_inferred;
info.trial_time = trial_time;
info.start_time = start_time;
info.choice_labels = choice_labels;
info.choice_frame = choice_frame;
info.trial_frame = trial_frame;
info.colorscan = colorscan;
info.mousedist = mousedist;
info.mouseID = mouseID;
if ~skipvid
info.bTeensyToCam = bTeensyToCam;
info.bCamToFrame = bCamToFrame;
end

% for python
info.sessionTime4Python = char(sessionTime);

end