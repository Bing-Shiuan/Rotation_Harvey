%% CALIBRATE COLOR

% This script reads the .dat behavior files and pulls out all the color
% info for each mouse.
% It is necessary to run before starting any social mice, and should be run
% every once in a while to make sure the thresholds are still calibrated

% ** modification of the script check_color_scan.mat

use_averaged = 0;
use_individual = 1;

addpath(genpath('../'))

%% Name of mice, in cell format, to check scans of
%miceNames = {'KM10','KM11','KM12','KM13'};
miceNames = {'KM10','KM11','KM13'};
%miceNames = {'KM16','KM18','KM19'}; %<- need to recolor before running today...
miceNames = {'KM20','KM21','KM22','KM23'}; % ugh the higher one isnt working that well...
% might need to try moving the lick port around, or lowering the scanner...
%miceNames = {'KM24','KM25','KM26','KM27'};
%miceNames = {'KM28','KM29','KM30','KM31'};
miceNames = {'KM33','KM34'};
miceNames = {'KM35','KM36'};
%miceNames = {'KM39','KM40','KM41','KM42'};

%miceNames = {'RF4','RF5','RF6','RF7'};

%miceNames = {'KM41','KM42','KM43','KM44'};

miceNames = {'KM45','KM46','KM47','KM48'};
% 
% 
% %miceNames = {'RF8','RF9','RF10'};
% 
% miceNames = {'KM49','KM50','KM51','KM52'};
%miceNames = {'KM53','KM54','KM55','KM56'};
%miceNames = {'KM57','KM58','KM59','KM60'};
%miceNames = {'KM61','KM62','KM63','KM64'};

miceNames = {'KM65','KM66','KM67'};
miceNames = {'KM49','KM50'};
%% Get color data for the most recent days

% labels
ports_labels = {'RT','RB','LT','LB'};
color_labels = {'Red','Green','Blue'};

portval_all = {};
color_data_all = [];
mouseIDall = [];
colormean = [];
color_unnorm = [];
choice_all = [];
mousedist_all = []; % this is what the teensy calculates.
portval_unnorm = {};

for m = 1:length(miceNames)
try
info = createBEHstruct_nonsocial(miceNames{m}, 'today', 'Color', 1);
catch
info = createBEHstruct_nonsocial(miceNames{m}, 2, 'Color', 1);
end

% if miceNames{m}=='KM67'
%     info = createBEHstruct_nonsocial(miceNames{m}, 3, 'Color', 1);
% end

choice = info.choice;
badid = cellfun(@isempty, choice);
colorscan = info.colorscan; 
% update badid for misscans???
badid = badid | sum(colorscan)<10 | sum(colorscan)>2000;
colorscan(:,badid) = []; choice(badid) = [];
choice = cellfun(@(v) v(1), choice);

color_unnorm = [color_unnorm, colorscan];
colornorm = colorscan(1:3, :)*255 ./ colorscan(4, :);
color_data_all = [color_data_all, colornorm];
mouseIDall = [mouseIDall, m*ones(1,length(colornorm))];
mousedist_all = [mousedist_all, info.mousedist];

portval = nan(3, 4);
portval_un = nan(4,4);
for j = unique(choice)
    portval(:,j+1) = mean(colorscan(1:3, choice==j)*255 ./ colorscan(4, choice==j), 2, 'omitnan');
    portval_un(:,j+1) = mean(colorscan(:, choice==j),2, 'omitnan');
end

portval_all{m} = portval;
colormean(:,m) = mean(colornorm,2,'omitnan');
choice_all = [choice_all, choice];
portval_unnorm{m} = portval_un;

end


ports_labels = info.choice_labels;


% determine threshold for what counts as a mouse?
% - based on pairwise distance between every single read?

% - calculate all pairwise color distnaces!
% pow(red - mouseRed,2) + pow(green - mouseGreen,2) + pow(blue-mouseBlue,2);

% a) caluclate within group distance
% b) calculate across group distance
% c) find the optimal threshold between the two!!c


% % better?
% figure; hold on; colorid = {'k','r','b','g'};ax1=[];
% for m = 1:4
% id = mouseIDall==m; coloruse = colorid{m};
% ax1(m) = scatter3(color_data_all(1,id), color_data_all(2,id), color_data_all(3,id), 5, coloruse);
% end
% legend(ax1, miceNames);
% xlabel('Red'); ylabel('Green'); zlabel('Blue');

% Compare distances to means
call = [color_data_all, colormean];
c = squareform(pdist(call', 'squaredeuclidean'));

if use_averaged
figure; hold on;
for jj = (length(miceNames)-1):-1:0
    plot(c(end-jj,:))
end
legend(miceNames)
xlabel('Trials'); ylabel('Squared color distance')
title('Color dist for global calibration');
end

% compare distances to mean portvals? and take the smallest?
% - so rather than compare to one mean, compare to multiple means...

% 
% % plot unnormalized scans?
% figure; hold on; colorid = {'k','r','b','g'};ax1=[];
% for m = 1:4
% id = mouseIDall==m; coloruse = colorid{m};
% ax1(m) = scatter3(color_unnorm(1,id), color_unnorm(2,id), color_unnorm(3,id), 5, coloruse);
% end
% legend(ax1, miceNames);
% xlabel('Red'); ylabel('Green'); zlabel('Blue');

% Instead of comparing to mean value, compare to individual port calibration
% replace nans with mean value in portval_all
portval_fill = portval_all;
for m = 1:length(portval_all)
    for portid = 1:4
        if any(isnan(portval_all{m}(:,portid)))
            portval_fill{m}(:,portid) = colormean(:,m);
        end
    end
end

colordist_by_port = zeros(length(portval_all),length(choice_all));
for m = 1:length(portval_all)
    for trial = 1:length(choice_all)
        colorsamp = color_data_all(:,trial);
        %colortest = portval_all{m}(:,choice_all(trial) + 1);
        colortest = portval_fill{m}(:,choice_all(trial) + 1);
        disttemp = pdist2(colorsamp', colortest', 'squaredeuclidean');
        colordist_by_port(m, trial) = disttemp;
    end
end

if use_individual
figure; hold on;
plot(colordist_by_port');
legend(miceNames)
xlabel('Trials'); ylabel('Squared color distance')
title('Color dist for individual calibrations');
end


if use_individual
% Print header
for m = 1:length(miceNames)
fprintf('\n')
fprintf(miceNames{m}); fprintf('\n');
fprintf('        ');
for c = 1:length(color_labels)
    fprintf('%10s', color_labels{c});
end
fprintf('\n');

% Print rows with data
data = round( portval_all{m}' );
for r = 1:length(ports_labels)
    fprintf('%-8s', ports_labels{r});  % Print row label
    for c = 1:size(data, 2)
        fprintf('%10.0f', data(r, c));  % Print data values
    end
    fprintf('\n');
end
fprintf('\n');
end
end

if use_averaged
data =  round( colormean' );
ports_labels = miceNames;
fprintf('\n');
fprintf('        ');
for c = 1:length(color_labels)
    fprintf('%10s', color_labels{c});
end
fprintf('\n');
for r = 1:length(ports_labels)
    fprintf('%-8s', ports_labels{r});  % Print row label
    for c = 1:size(data, 2)
        fprintf('%10.0f', data(r, c));  % Print data values
    end
    fprintf('\n');
end
fprintf('\n');
end

%% scan many days back to get the port value colors?
if 0

portval_day = {}; numchoice_day = {};
m = 3; % which mouse to use

for days_back = 0:30;
try
info = createBEHstruct_nonsocial(miceNames{m}, days_back, 'Color', 1);
catch
    %disp(days_back)
    continue;
end
choice = info.choice;
badid = cellfun(@isempty, choice);
colorscan = info.colorscan; 
% update badid for misscans???
badid = badid | sum(colorscan)<10 | sum(colorscan)>2000;
colorscan(:,badid) = []; choice(badid) = [];
choice = cellfun(@(v) v(1), choice);

portval = nan(3, 4);
for j = unique(choice)
    portval(:,j+1) = mean(colorscan(1:3, choice==j)*255 ./ colorscan(4, choice==j), 2);
end

portval_day{end+1} = portval;
numchoice_day{end+1} = choice;

end

end