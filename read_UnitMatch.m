clc;clear;close all
matching_mouse=["KM62" "KM61-62" "KM62" "KM61-62"]
matching_dates=["251108" "251109" "251112" "251113"]
unitMatchPath="E:\BingShiuan\KM62Bottom"
UM_Prob_th=0.9
%%
MatchTable=readtable(strcat(unitMatchPath,"\MatchTable.csv"));
%%
matched_pairs=MatchTable(MatchTable.UMProbabilities>UM_Prob_th ...
    & MatchTable.ID1~=MatchTable.ID2 ...
    & MatchTable.RecSes1~=MatchTable.RecSes2,:);
matchingMap=repmat([-1 -1 -1 -1],height(matched_pairs),1);
for p=1:height(matched_pairs)
    matchingMap(p,matched_pairs.RecSes1(p))=matched_pairs.ID1(p);
    matchingMap(p,matched_pairs.RecSes2(p))=matched_pairs.ID2(p);
end

        
