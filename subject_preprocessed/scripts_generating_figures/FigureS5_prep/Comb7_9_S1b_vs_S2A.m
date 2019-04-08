tmQ=[0;6;9;12;15;22];

%% Perc = s1b%
Subset_7 = Quantify(Quantify(:,6)==7,:);
figure(1); hold on;
Ratio_7 = (Subset_7(:,3)+Subset_7(:,4))./(Subset_7(:,1)+Subset_7(:,2)+Subset_7(:,3)+Subset_7(:,4));
Ratio_7_mean_ts = []; Ratio_7_std_ts = [];
for i = 0:5;
    Ratio_7_mean_ts = [Ratio_7_mean_ts; mean(Ratio_7(find(Subset_7(:,5)==i)))]; 
    Ratio_7_std_ts = [Ratio_7_std_ts; var(Ratio_7(find(Subset_7(:,5)==i)))];
end
errorbar(tmQ, Ratio_7_mean_ts,Ratio_7_std_ts,'r');
ylim([0 1])
%%
Subset_8 = Quantify(Quantify(:,6)==8,:);
Ratio_8 = (Subset_8(:,3)+Subset_8(:,4))./(Subset_8(:,1)+Subset_8(:,2)+Subset_8(:,3)+Subset_8(:,4));
Ratio_8_mean_ts = []; Ratio_8_std_ts = [];
for i = 0:5;
    Ratio_8_mean_ts = [Ratio_8_mean_ts; mean(Ratio_8(find(Subset_8(:,5)==i)))]; 
    Ratio_8_std_ts = [Ratio_8_std_ts; var(Ratio_8(find(Subset_8(:,5)==i)))];
end
errorbar(tmQ, Ratio_8_mean_ts,Ratio_8_std_ts,'r');
ylim([0 1])
    
%%
Subset_9 = Quantify(Quantify(:,6)==9,:);
Ratio_9 = (Subset_9(:,3)+Subset_9(:,4))./(Subset_9(:,1)+Subset_9(:,2)+Subset_9(:,3)+Subset_9(:,4));
Ratio_9_mean_ts = []; Ratio_9_std_ts = [];
for i = 0:5;
    Ratio_9_mean_ts = [Ratio_9_mean_ts; mean(Ratio_9(find(Subset_9(:,5)==i)))]; 
    Ratio_9_std_ts = [Ratio_9_std_ts; var(Ratio_9(find(Subset_9(:,5)==i)))];
end
errorbar(tmQ, Ratio_9_mean_ts,Ratio_9_std_ts,'r');
ylim([0 1])
    
    