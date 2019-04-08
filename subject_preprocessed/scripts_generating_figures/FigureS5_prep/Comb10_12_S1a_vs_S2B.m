tmQ=[0;6;9;12;15;22];

%% Perc = s1b%
Subset_10 = Quantify(Quantify(:,6)==10,:);
figure(1); hold on;
Ratio_10 = (Subset_10(:,3)+Subset_10(:,4))./(Subset_10(:,1)+Subset_10(:,2)+Subset_10(:,3)+Subset_10(:,4));
Ratio_10_mean_ts = []; Ratio_10_std_ts = [];
for i = 0:5;
    Ratio_10_mean_ts = [Ratio_10_mean_ts; mean(Ratio_10(find(Subset_10(:,5)==i)))]; 
    Ratio_10_std_ts = [Ratio_10_std_ts; var(Ratio_10(find(Subset_10(:,5)==i)))];
end
errorbar(tmQ, Ratio_10_mean_ts,Ratio_10_std_ts,'b');
ylim([0 1])
%%
Subset_11 = Quantify(Quantify(:,6)==11,:);
Ratio_11 = (Subset_11(:,3)+Subset_11(:,4))./(Subset_11(:,1)+Subset_11(:,2)+Subset_11(:,3)+Subset_11(:,4));
Ratio_11_mean_ts = []; Ratio_11_std_ts = [];
for i = 0:5;
    Ratio_11_mean_ts = [Ratio_11_mean_ts; mean(Ratio_11(find(Subset_11(:,5)==i)))]; 
    Ratio_11_std_ts = [Ratio_11_std_ts; var(Ratio_11(find(Subset_11(:,5)==i)))];
end
errorbar(tmQ, Ratio_11_mean_ts,Ratio_11_std_ts,'b');
ylim([0 1])
    
%%
Subset_12 = Quantify(Quantify(:,6)==12,:);
Ratio_12 = (Subset_12(:,3)+Subset_12(:,4))./(Subset_12(:,1)+Subset_12(:,2)+Subset_12(:,3)+Subset_12(:,4));
Ratio_12_mean_ts = []; Ratio_12_std_ts = [];
for i = 0:5;
    Ratio_12_mean_ts = [Ratio_12_mean_ts; mean(Ratio_12(find(Subset_12(:,5)==i)))]; 
    Ratio_12_std_ts = [Ratio_12_std_ts; var(Ratio_12(find(Subset_12(:,5)==i)))];
end
errorbar(tmQ, Ratio_12_mean_ts,Ratio_12_std_ts,'b');
ylim([0 1])
    
    