tmQ=[0;6;9;12;15;22];

%% Perc = s1b%
Subset_13 = Quantify(Quantify(:,6)==13,:);
figure(1); hold on;
Ratio_13 = (Subset_13(:,3)+Subset_13(:,4))./(Subset_13(:,1)+Subset_13(:,2)+Subset_13(:,3)+Subset_13(:,4));
Ratio_13_mean_ts = []; Ratio_13_std_ts = [];
for i = 0:5;
    Ratio_13_mean_ts = [Ratio_13_mean_ts; mean(Ratio_13(find(Subset_13(:,5)==i)))]; 
    Ratio_13_std_ts = [Ratio_13_std_ts; var(Ratio_13(find(Subset_13(:,5)==i)))];
end
errorbar(tmQ, Ratio_13_mean_ts,Ratio_13_std_ts,'r');
ylim([0 1])
%%
Subset_14 = Quantify(Quantify(:,6)==14,:);
Ratio_14 = (Subset_14(:,3)+Subset_14(:,4))./(Subset_14(:,1)+Subset_14(:,2)+Subset_14(:,3)+Subset_14(:,4));
Ratio_14_mean_ts = []; Ratio_14_std_ts = [];
for i = 0:5;
    Ratio_14_mean_ts = [Ratio_14_mean_ts; mean(Ratio_14(find(Subset_14(:,5)==i)))]; 
    Ratio_14_std_ts = [Ratio_14_std_ts; var(Ratio_14(find(Subset_14(:,5)==i)))];
end
errorbar(tmQ, Ratio_14_mean_ts,Ratio_14_std_ts,'r');
ylim([0 1])
    
%%
Subset_15 = Quantify(Quantify(:,6)==15,:);
Ratio_15 = (Subset_15(:,3)+Subset_15(:,4))./(Subset_15(:,1)+Subset_15(:,2)+Subset_15(:,3)+Subset_15(:,4));
Ratio_15_mean_ts = []; Ratio_15_std_ts = [];
for i = 0:5;
    Ratio_15_mean_ts = [Ratio_15_mean_ts; mean(Ratio_15(find(Subset_15(:,5)==i)))]; 
    Ratio_15_std_ts = [Ratio_15_std_ts; var(Ratio_15(find(Subset_15(:,5)==i)))];
end
errorbar(tmQ, Ratio_15_mean_ts,Ratio_15_std_ts,'r');
ylim([0 1])
    
    