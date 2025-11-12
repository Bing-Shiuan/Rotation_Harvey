% clc;clear;close all;
mouseName="KM49-50";
date="251104";

pathCSV=dir(strcat("Z:\HarveyLab\Tier1\Kevin\Videos\",mouseName,"\",date,"*snapshot_best-530.csv"));
inputCSV=strcat(pathCSV.folder,"\",pathCSV.name)
opts = detectImportOptions(inputCSV,'NumHeaderLines',2);
tic
T = readtable(inputCSV,opts);
disp("CSV loaded")
toc
XX=["x" "x_1"];
YY=["y" "y_1"];

TT=T;

tic
for i=1:length(XX)
    disp(strcat("filtering animal",string(i),"..."))
v =[0;sqrt(diff(T.(XX(i))).^2+diff(T.(YY(i))).^2)];
filtered_x= smooth(T.(XX(i)),10,'rloess');
filtered_y= smooth(T.(YY(i)),10,'rloess');
filtered_v =[0;sqrt(diff(filtered_x).^2+diff(filtered_y).^2)];
figure;
subplot(3,1,1)
plot(T.(XX(i)))
hold on
plot(filtered_x)
xlim([1 2500])
title("x")


subplot(3,1,2)
plot(T.(YY(i)))
hold on
plot(filtered_y)
xlim([1 2500])
title("y")

subplot(3,1,3)
plot(v)
hold on
plot(filtered_v)
xlim([1 2500])
title("v")

TT.(XX(i))=filtered_x;
TT.(YY(i))=filtered_y;
end
toc
% %% save coordinates
save(strcat("Z:\HarveyLab\Tier1\Bing_Shiuan\Codes\",mouseName,"_",date,"_tracking.mat"),"TT")
%%
% --- minimal labeled video writer (no functions) ---
pathVideo = dir(strcat("Z:\HarveyLab\Tier1\Kevin\Videos\",mouseName,"\",date,"*.mp4"));
inputVideo  = strcat(pathVideo.folder,"\",pathVideo.name)
outputVideo = strcat(pathVideo.folder,"\refined_",pathVideo.name)
xA=TT.x;
yA=TT.y;
xB=TT.x_1;
yB=TT.y_1;
% Arrays must exist: xA, yA, xB, yB (length = #frames). Use NaN to skip a point.
% If your coordinates are normalized [0,1], set this true to scale to pixels:
useNormalized = false;

v = VideoReader(inputVideo);
w = VideoWriter(outputVideo, 'MPEG-4');
w.FrameRate = v.FrameRate;
open(w);
% Total frames based on arrays (safer than Duration*FPS for some codecs)
N = min([numel(xA), numel(yA), numel(xB), numel(yB)]);
h = waitbar(0, sprintf('Processing 0/%d frames...', N), 'Name','Annotating video');
cleanupWait = onCleanup(@() (ishandle(h) && delete(h)));  %#ok<NASGU>
f = 0;
while hasFrame(v)
    I = readFrame(v);
    f = f + 1;

    if f > numel(xA) || f > numel(xB)
        % Safety stop if arrays are shorter than video (shouldn't happen per your note)
        break;
    end

    [H,W,~] = size(I);

    % Mouse A
    xa = xA(f); ya = yA(f);
    if useNormalized
        xa = xa * W; ya = ya * H;
    end
    if ~isnan(xa) && ~isnan(ya)
        xa = min(max(xa,1), W); ya = min(max(ya,1), H);
        I = insertMarker(I, [xa ya], 'o', 'Color', 'yellow', 'Size', 10);
        I = insertText(I, [xa+8 ya-18], 'A', 'FontSize', 16, ...
                       'BoxOpacity', 0.6, 'TextColor', 'white', 'BoxColor', 'yellow');
    end

    % Mouse B
    xb = xB(f); yb = yB(f);
    if useNormalized
        xb = xb * W; yb = yb * H;
    end
    if ~isnan(xb) && ~isnan(yb)
        xb = min(max(xb,1), W); yb = min(max(yb,1), H);
        I = insertMarker(I, [xb yb], 'o', 'Color', 'cyan', 'Size', 10);
        I = insertText(I, [xb+8 yb-18], 'B', 'FontSize', 16, ...
                       'BoxOpacity', 0.6, 'TextColor', 'white', 'BoxColor', 'cyan');
    end

    writeVideo(w, I);

        % ---- progress bar (update every ~10 frames) ----
    if ishandle(h) && (mod(f,10)==0 || f==N)
        frac = min(1, f/N);
        waitbar(frac, h, sprintf('Processing %d/%d (%.0f%%)', f, N, 100*frac));
        drawnow;  % ensure UI updates
    end
end

close(w);
disp(['✅ Wrote ', outputVideo]);
