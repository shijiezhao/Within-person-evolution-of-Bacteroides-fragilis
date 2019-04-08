tmQ=[0;6;9;12;15;22];

%% Perc = s1b%
Subset_19_1 = Quantify(Quantify(:,6)==19 & Quantify(:,7)==1,:);
Subset_19_2 = Quantify(Quantify(:,6)==19 & Quantify(:,7)==2,:);

figure(1); hold on;
% This is the perc of s1
Ratio_19_1 = (Subset_19_1(:,3)+Subset_19_1(:,4))./(Subset_19_1(:,1)+Subset_19_1(:,2)+Subset_19_1(:,3)+Subset_19_1(:,4));
% perc of s1-a-1-1
Ratio_19_2 = (Subset_19_2(:,1)+Subset_19_2(:,2))./(Subset_19_2(:,1)+Subset_19_2(:,2)+Subset_19_2(:,3)+Subset_19_2(:,4));

Ratio_19_1_mean_ts = []; Ratio_19_1_std_ts = [];
for i = 0:5;
    Ratio_19_1_mean_ts = [Ratio_19_1_mean_ts; mean(Ratio_19_1(find(Subset_19_1(:,5)==i)))]; 
    Ratio_19_1_std_ts = [Ratio_19_1_std_ts; var(Ratio_19_1(find(Subset_19_1(:,5)==i)))];
end
errorbar(tmQ, 1-Ratio_19_1_mean_ts,Ratio_19_1_std_ts,'r');
ylim([0 1])

Ratio_19_2_mean_ts = []; Ratio_19_2_std_ts = [];
for i = 0:5;
    Ratio_19_2_mean_ts = [Ratio_19_2_mean_ts; mean(Ratio_19_2(find(Subset_19_2(:,5)==i)))]; 
    Ratio_19_2_std_ts = [Ratio_19_2_std_ts; var(Ratio_19_2(find(Subset_19_2(:,5)==i)))];
end
errorbar(tmQ, 1-Ratio_19_2_mean_ts,Ratio_19_2_std_ts,'b');
ylim([0 1])

Ratio_19_3_mean_ts = Ratio_19_1_mean_ts - Ratio_19_2_mean_ts;
errorbar(tmQ, Ratio_19_3_mean_ts,Ratio_19_2_std_ts,'k');
ylim([0 1])

errorbar(tmQ, Ratio_19_3_mean_ts./(1-Ratio_19_2_mean_ts),Ratio_19_2_std_ts,'g');
ylim([0 1])
