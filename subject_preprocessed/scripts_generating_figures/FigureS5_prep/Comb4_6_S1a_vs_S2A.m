tmQ=[0;6;9;12;15;22];

%% Perc = s1b%
Subset_4 = Quantify(Quantify(:,6)==4,:);
figure(1); hold on;
Ratio_4 = (Subset_4(:,3)+Subset_4(:,4))./(Subset_4(:,1)+Subset_4(:,2)+Subset_4(:,3)+Subset_4(:,4));
Ratio_4_mean_ts = []; Ratio_4_std_ts = [];
for i = 0:5;
    Ratio_4_mean_ts = [Ratio_4_mean_ts; mean(Ratio_4(find(Subset_4(:,5)==i)))]; 
    Ratio_4_std_ts = [Ratio_4_std_ts; var(Ratio_4(find(Subset_4(:,5)==i)))];
end
errorbar(tmQ, Ratio_4_mean_ts,Ratio_4_std_ts,'b');
Ratios{4} = Ratio_4_mean_ts;
ylim([0 1])
%%
Subset_5 = Quantify(Quantify(:,6)==5,:);
Ratio_5 = (Subset_5(:,3)+Subset_5(:,4))./(Subset_5(:,1)+Subset_5(:,2)+Subset_5(:,3)+Subset_5(:,4));
Ratio_5_mean_ts = []; Ratio_5_std_ts = [];
for i = 0:5;
    Ratio_5_mean_ts = [Ratio_5_mean_ts; mean(Ratio_5(find(Subset_5(:,5)==i)))]; 
    Ratio_5_std_ts = [Ratio_5_std_ts; var(Ratio_5(find(Subset_5(:,5)==i)))];
end
errorbar(tmQ, Ratio_5_mean_ts,Ratio_5_std_ts,'b');
Ratios{5} = Ratio_5_mean_ts;
ylim([0 1])
    
%%
Subset_6 = Quantify(Quantify(:,6)==6,:);
Ratio_6 = (Subset_6(:,3)+Subset_6(:,4))./(Subset_6(:,1)+Subset_6(:,2)+Subset_6(:,3)+Subset_6(:,4));
Ratio_6_mean_ts = []; Ratio_6_std_ts = [];
for i = 0:5;
    Ratio_6_mean_ts = [Ratio_6_mean_ts; mean(Ratio_6(find(Subset_6(:,5)==i)))]; 
    Ratio_6_std_ts = [Ratio_6_std_ts; var(Ratio_6(find(Subset_6(:,5)==i)))];
end
errorbar(tmQ, Ratio_6_mean_ts,Ratio_6_std_ts,'b');
Ratios{6} = Ratio_6_mean_ts;

ylim([0 1])
    
    