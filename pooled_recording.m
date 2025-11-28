%% mean Fr cross days
% clc; clear; close all;
date= '251109'
mouse_name = 'KM61-62' %**0 1 swap
load(strcat("Z:\HarveyLab\Tier1\Bing_Shiuan\Codes\",mouse_name,"_",date,"aligned.mat"))
FR_all_3Darray=[];
for i=1:length(unit_depth)
FR_all_3Darray(:,:,i)=FR_all.(strcat("unit",string(i)));
end
baseFR=permute(mean(FR_all_3Darray(1:30,:,:),1),[3 2 1]);
goodunit=find(goodid)
%%
figure;
plot(baseFR(245,:))
%%
zs_baseFr=zscore(baseFR,0,2);
[~,I]=sort(mean(zs_baseFr(:,1:30)-zs_baseFr(:,end-29:end),2),"descend");

figure;
subplot(3,1,1)
heatmap(zs_baseFr(I,:),"GridVisible","off")
clim([-2 4])
colormap("turbo")
title(date)
subplot(3,1,2)
plot(mean(baseFR,1))
hold on
yyaxis right
plot(ITI)
ylim([1 600000])
subplot(3,1,3)
plot(cumsum(ITI))
%%
% figure;
% plot(mean(baseFR,1))%all_unit_MeanFR
% hold on
% yyaxis right
% plot(mean(abs(trial_mean_amp),1))
% Q, dQ, ddQ
deltaQ=Q1(1,:)-Q1(2,:);
deltaQ2=Q2(1,:)-Q2(2,:);
diff_deltaQ = [0;diff(smooth(deltaQ,20))]';
diff_deltaQ2 = [0;diff(smooth(deltaQ2,20))]';
figure;
plot(deltaQ);hold on
plot(deltaQ2)
Corr_Q=corr(deltaQ',deltaQ2');
dissociate_trials=find(deltaQ+deltaQ2>-2);
Fr_use=baseFR(:,dissociate_trials);

figure;

trial_param=deltaQ(:,dissociate_trials);
Corr_to_param=corr(Fr_use',trial_param');
[~,Iunit]=sort(Corr_to_param,1,"descend");
[Strial,IdQ]=sort(trial_param,2,"descend");

subplot(3,2,1)
plot(Strial)
xlim([1 length(Strial)])


baseFR_sort_trial=Fr_use(:,IdQ);
baseFR_sort_both=baseFR_sort_trial(Iunit,:);

subplot(3,2,3)
heatmap(zscore(baseFR_sort_both,0,2),"GridVisible","off")
colormap("turbo")
clim([-0.5 2.5])
subplot(3,2,5)
histogram(Corr_to_param,[-0.75:0.1:0.75],"FaceAlpha",0.5,"Normalization","percentage")
hold on


trial_param=deltaQ2(:,dissociate_trials);
Corr_to_param=corr(Fr_use',trial_param');

[~,Iunit]=sort(Corr_to_param,1,"descend");
[Strial,IdQ]=sort(trial_param,2,"descend");
baseFR_sort_trial=Fr_use(:,IdQ);
baseFR_sort_both=baseFR_sort_trial(Iunit,:);

subplot(3,2,2)
plot(Strial)
xlim([1 length(Strial)])

subplot(3,2,4)
heatmap(zscore(baseFR_sort_both,0,2),"GridVisible","off")
colormap("turbo")
clim([-0.5 2.5])

subplot(3,2,5)
histogram(Corr_to_param,[-0.75:0.1:0.75],"FaceAlpha",0.5,"Normalization","percentage")

subplot(3,2,6)
scatter(deltaQ,deltaQ2)
hold on
scatter(deltaQ(dissociate_trials),deltaQ2(dissociate_trials))
title(string(Corr_Q))
%% 4 condition mean
% clc; clear; close all;
% date= ["251104" "251106"];
date= ["251021" "251023"];
mouse_name = 'KM49-50' %**0 1 swap
trial_avg_all_days=[];
 unit_feature_all_days=[];
unit_number_record=[];
unit_medlat_all=[];
unit_depth_all=[];
for d=date
load(strcat("Z:\HarveyLab\Tier1\Bing_Shiuan\Codes\",mouse_name,"_",d,"aligned.mat"))
zs_FR_3Darray=[];
unit_medlat_all=[unit_medlat_all;unit_medlat];
unit_depth_all=[unit_depth_all;unit_depth];
for i=1:length(unit_depth)
zs_FR_3Darray(:,:,i)=zscore(FR_all.(strcat("unit",string(i))),0,1);
end
zs_FR_3Darray=permute(zs_FR_3Darray,[1 3 2]);


%
condition_c1=[-1 1];
condition_c2=[-1 1];
condition_r1=[-1 1];
condition_r2=[-1 1];
unit_feature_of_cond=[];
unit_feature_all=[];
trial_avg_all=[];
Condition_label_all=[];
 % for c1=condition_c1
 %     for c2=condition_c2
 %         trial_ID=find(choice1==c1 & choice2==c2);
 %          Condition_label=strcat("c1:",string(c1), " c2:",string(c2));
 for r1=condition_r1 %condition_r1
     for  r2=condition_r2%condition_r2
         trial_ID=find(reward1==r1 & reward2==r2);
         Condition_label=strcat("r1:",string(r1), " r2:",string(r2));
         trial_avg=mean(zs_FR_3Darray(:,:,trial_ID),3);
         unit_feature_of_cond=[mean(trial_avg(21:50,:),1);mean(trial_avg(51:65,:),1);mean(trial_avg(66:95,:),1)];
         trial_avg_all=cat(1,trial_avg_all,trial_avg);
         unit_feature_all=[unit_feature_all;unit_feature_of_cond];
         Condition_label_all=[Condition_label_all;[strcat(Condition_label,":-2 to 0");strcat(Condition_label,":0 to 1");strcat(Condition_label,":1 to 3")]];

     end
 end
 trial_avg_all_days=cat(2,trial_avg_all_days,trial_avg_all);%concate in unit
 unit_feature_all_days=cat(2,unit_feature_all_days,unit_feature_all);
 unit_number_record=[unit_number_record length(unit_depth)];
end
 [coeff_cl,score_cl,~,~,explained_cl,~] = pca(unit_feature_all_days');

 figure;
 bar(explained_cl)
 hold on
 plot(cumsum(explained_cl),'-o')

%% cluster by...
numClusters=4;
clustering_method=score_cl(:,1:3);% here
linkage_method="ward";
%clustering and heatmap plotting
Z = linkage(clustering_method,linkage_method);
cg_align = clustergram(clustering_method,'Colormap',turbo,'Cluster', 'column'); %clustering method
% cg_align.DisplayRange =8
set(cg_align,'Linkage',linkage_method,'Dendrogram',12)
indx= str2double(string(cg_align.RowLabels)); %take the index order of clustering because clustergram uses "optimal leaf order" that maximize the similarity of adjacent leaves
T = cluster(Z,'maxclust',numClusters,'Criterion','distance');

%% cluster waveform
figure;
for cl=1:max(T)
    subplot(max(T),1,cl)
    plot([-4:0.0665:6],mean(trial_avg_all_days(1:151,find(T==cl)),2));
    hold on
    plot([-4:0.0665:6],mean(trial_avg_all_days(152:302,find(T==cl)),2))
    hold on
    plot([-4:0.0665:6],mean(trial_avg_all_days(303:453,find(T==cl)),2))
    hold on
    plot([-4:0.0665:6],mean(trial_avg_all_days(454:604,find(T==cl)),2))
    title(strcat("n=",string(length(find(T==cl)))))
    
end
%% position
unit_number=[0 unit_number_record];
for i=1:length(unit_number_record)
 figure;
 scatter(unit_medlat_all(unit_number(i)+1:unit_number(i+1)+unit_number(i)),unit_depth_all(unit_number(i)+1:unit_number(i+1)+unit_number(i)), ...
     30,T(unit_number(i)+1:unit_number(i+1)+unit_number(i)),"filled","XJitter","density",XJitterWidth = 80)
colormap("turbo")
clim([0.5 max(T)+0.5])
 xlim([-100 800])
 title(date(i))
end

 %%
 figure;
 scatter3(score_cl(:,1),score_cl(:,2),score_cl(:,3),[],T,"filled")
 xlabel("PC1")
 ylabel("PC2")
 zlabel("PC3")
colormap("turbo")
clim([0.5 max(T)+0.5])


%% cluster

%% position
