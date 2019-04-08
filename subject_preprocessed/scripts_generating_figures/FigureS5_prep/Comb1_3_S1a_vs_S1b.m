%% Timepoint for figure S5G-S5H
tmQ=[0;6;9;12;15;22];

%% About Quantify:
% First and Third columns: reads aligned to the ancestral allele
% Second and Fourth columns: reads aligned to the mutant allele
% Fifth column: time point
% Sixth column: combinations
% Seventh column: different types of amplicon
% Eighth column: technical replicates
% Combinations:
% #1-#3: 

%% Perc = s1b%
Subset_1 = Quantify(Quantify(:,6)==1,:);
figure(1); hold on;
Ratio_1 = (Subset_1(:,1)+Subset_1(:,2))./(Subset_1(:,1)+Subset_1(:,2)+Subset_1(:,3)+Subset_1(:,4));
Ratio_1_mean_ts = []; Ratio_1_std_ts = [];
for i = 0:5;
    Ratio_1_mean_ts = [Ratio_1_mean_ts; mean(Ratio_1(find(Subset_1(:,5)==i)))]; 
    Ratio_1_std_ts = [Ratio_1_std_ts; var(Ratio_1(find(Subset_1(:,5)==i)))];
end
errorbar(tmQ, Ratio_1_mean_ts,Ratio_1_std_ts);

ylim([0 1])
%%
Subset_2 = Quantify(Quantify(:,6)==2,:);
Ratio_2 = (Subset_2(:,1)+Subset_2(:,2))./(Subset_2(:,1)+Subset_2(:,2)+Subset_2(:,3)+Subset_2(:,4));
Ratio_2_mean_ts = []; Ratio_2_std_ts = [];
for i = 0:5;
    Ratio_2_mean_ts = [Ratio_2_mean_ts; mean(Ratio_2(find(Subset_2(:,5)==i)))]; 
    Ratio_2_std_ts = [Ratio_2_std_ts; var(Ratio_2(find(Subset_2(:,5)==i)))];
end
errorbar(tmQ, Ratio_2_mean_ts,Ratio_2_std_ts);
ylim([0 1])
    
%%
Subset_3 = Quantify(Quantify(:,6)==3,:);
Ratio_3 = (Subset_3(:,1)+Subset_3(:,2))./(Subset_3(:,1)+Subset_3(:,2)+Subset_3(:,3)+Subset_3(:,4));
Ratio_3_mean_ts = []; Ratio_3_std_ts = [];
for i = 0:5;
    Ratio_3_mean_ts = [Ratio_3_mean_ts; mean(Ratio_3(find(Subset_3(:,5)==i)))]; 
    Ratio_3_std_ts = [Ratio_3_std_ts; var(Ratio_3(find(Subset_3(:,5)==i)))];
end
errorbar(tmQ, Ratio_3_mean_ts,Ratio_3_std_ts);
ylim([0 1])
    
 