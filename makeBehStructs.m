%% Loop through mice to make info sessions for analyses

addpath(genpath('../'));

data_path = 'Z:\HarveyLab\Tier1\Kevin\Videos\'; % <- this points to where data is saved



% Name of mice to save
miceNonsocial = {'KM33','KM34','KM35','KM36','KM37','KM38'};

miceSocial= {'KM33-34','KM35-36','KM37-38'};
% 
miceNonsocial = {'KM41','KM42'};
miceSocial = { 'KM41-42'};
 
miceNonsocial = {'KM43','KM44'};

miceNonsocial = {'KM41','KM42','KM43','KM44','KM45','KM46','KM47','KM48'};
miceSocial = {'KM41-42','KM43-44','KM45-46','KM47-48'};


miceNonsocial = {'KM49','KM50','KM51','KM52','KM53','KM54','KM55','KM56'};
miceSocial = {};

miceNonsocial = {'KM49','KM50'};
miceSocial = {'KM49-50'};

boxtype = 'T-maze2';


%% 
% data_path = 'Z:\HarveyLab\Tier1\Rhyanne\';
% miceNonsocial = {'RF1','RF2','RF3','RF4','RF5','RF6','RF7'};
% miceSocial = {'RF4-6','RF5-7'};
% 
% miceNonsocial = {'RF1','RF2','RF3','RF4','RF5','RF6','RF7','RF8','RF9','RF10','RF11','RF12','RF13','RF14'};
% miceSocial = {'RF1-2','RF1-3','RF2-3','RF4-6','RF5-7','RF8-9','RF9-10','RF8-10','RF11-12','RF13-14'};
% miceNonsocial = {'RF8','RF9','RF10'};
% miceSocial = {'RF8-9','RF9-10','RF8-10'};
% 
% boxtype = 'T-maze2';

%% Non social

for m = 1:length(miceNonsocial)

    mouse_name = miceNonsocial{m};
    disp(mouse_name)
    mouse_files = dir(fullfile(data_path, mouse_name, '*.dat'));
    mouse_files([mouse_files.bytes]<100000) = []; %<- need to filter better!
    % check if multipel files in a day, and remove the smallest?
    fnames = {}; for s = 1:length(mouse_files); fnames{s} = mouse_files(s).name(1:6); end
    badid = [];
    for f = 1:length(fnames)
        if sum(strcmp(fnames, fnames{f}))>1
            id = find(strcmp(fnames, fnames{f}));
            [~,keepid] = max([mouse_files(id).bytes]);
            id(keepid) = []; 
            badid = [badid, id];
        end
    end
    badid = unique(badid);
    mouse_files(badid) = [];

    clear info
    
    for sess = 1:length(mouse_files)
        mouse_path = fullfile(mouse_files(sess).folder, mouse_files(sess).name);
        info(sess) = createBEHstruct_nonsocial(mouse_name, mouse_path, boxtype);
    end
    savepath = fullfile(data_path, mouse_name, 'mouseBEHstruct.mat');
    save(savepath, 'info')

end


%% Social



for m = 1:length(miceSocial)
    disp(m)
    mouse_name = miceSocial{m};
    mouse_files = dir(fullfile(data_path, mouse_name, '*.dat'));
    mouse_files([mouse_files.bytes]<100000) = []; %<- need to filter better!
    % check if multipel files in a day, and remove the smallest?
    fnames = {}; for s = 1:length(mouse_files); fnames{s} = mouse_files(s).name(1:6); end
    badid = [];
    for f = 1:length(fnames)
        if sum(strcmp(fnames, fnames{f}))>1
            % if multiple files in one day, take only largest file (fake
            % starts)
            id = find(strcmp(fnames, fnames{f}));
            [~,keepid] = max([mouse_files(id).bytes]);
            id(keepid) = []; 
            badid = [badid, id];
        end
    end
    badid = unique(badid);
    mouse_files(badid) = [];

    clear info
    
    for sess = 1:length(mouse_files)
        mouse_path = fullfile(mouse_files(sess).folder, mouse_files(sess).name);
        info(sess) = createBEHstruct_social(mouse_name, mouse_path, boxtype);
    end
    savepath = fullfile(data_path, mouse_name, 'mouseBEHstruct.mat');
    save(savepath, 'info')

end


