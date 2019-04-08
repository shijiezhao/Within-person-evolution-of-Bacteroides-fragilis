Ratio_16_4 = Ratio_16_2./Ratio_16_1;
Ratio_16_4_mean_ts = []; Ratio_16_4_std_ts = [];
clear std
for i = 0:5;
    tmp = [];
    for j = 1:3;
        if i ==0; jj=1;else;jj=j;end
        TMP1 = Ratio_16_1(find(Subset_16_1(:,5)==i & Subset_16_1(:,8)==jj));
        TMP2 = Ratio_16_2(find(Subset_16_2(:,5)==i & Subset_16_2(:,8)==jj));
        tmp = [tmp;(TMP1-TMP2)/(1-TMP2)];
    end
    Ratio_16_4_mean_ts = [Ratio_16_4_mean_ts;mean(tmp)]; 
    Ratio_16_4_std_ts = [Ratio_16_4_std_ts;std(tmp)];
end
%%
Ratio_17_4 = Ratio_17_2./Ratio_17_1;
Ratio_17_4_mean_ts = []; Ratio_17_4_std_ts = [];

for i = 0:5;
    tmp = [];
    for j = 1:3;
        if i ==0; jj=1;else;jj=j;end
        TMP1 = Ratio_17_1(find(Subset_17_1(:,5)==i & Subset_17_1(:,8)==jj));
        TMP2 = Ratio_17_2(find(Subset_17_2(:,5)==i & Subset_17_2(:,8)==jj));
        tmp = [tmp;(TMP1-TMP2)/(1-TMP2)];
    end
    Ratio_17_4_mean_ts = [Ratio_17_4_mean_ts;mean(tmp)]; 
    Ratio_17_4_std_ts = [Ratio_17_4_std_ts;std(tmp)];
end
%%
Ratio_18_4 = Ratio_18_2./Ratio_18_1;
Ratio_18_4_mean_ts = []; Ratio_18_4_std_ts = [];

for i = 0:5;
    tmp = [];
    for j = 1:3;
        if i ==0; jj=1;else;jj=j;end
        TMP1 = Ratio_18_1(find(Subset_18_1(:,5)==i & Subset_18_1(:,8)==jj));
        TMP2 = Ratio_18_2(find(Subset_18_2(:,5)==i & Subset_18_2(:,8)==jj));
        tmp = [tmp;(TMP1-TMP2)/(1-TMP2)];
    end
    Ratio_18_4_mean_ts = [Ratio_18_4_mean_ts;mean(tmp)]; 
    Ratio_18_4_std_ts = [Ratio_18_4_std_ts;std(tmp)];
end
%%
Ratio_19_4 = Ratio_19_2./Ratio_19_1;
Ratio_19_4_mean_ts = []; Ratio_19_4_std_ts = [];

for i = 0:5;
    tmp = [];
    for j = 1:3;
        if i ==0; jj=1;else;jj=j;end
        TMP1 = Ratio_19_1(find(Subset_19_1(:,5)==i & Subset_19_1(:,8)==jj));
        TMP2 = Ratio_19_2(find(Subset_19_2(:,5)==i & Subset_19_2(:,8)==jj));
        tmp = [tmp;(TMP1-TMP2)/(1-TMP2)];
    end
    Ratio_19_4_mean_ts = [Ratio_19_4_mean_ts;mean(tmp)]; 
    Ratio_19_4_std_ts = [Ratio_19_4_std_ts;std(tmp)];
end
%%
Ratio_20_4 = Ratio_20_2./Ratio_20_1;
Ratio_20_4_mean_ts = []; Ratio_20_4_std_ts = [];

for i = 0:5;
    tmp = [];
    for j = 1:3;
        if i ==0; jj=1;else;jj=j;end
        TMP1 = Ratio_20_1(find(Subset_20_1(:,5)==i & Subset_20_1(:,8)==jj));
        TMP2 = Ratio_20_2(find(Subset_20_2(:,5)==i & Subset_20_2(:,8)==jj));
        tmp = [tmp;(TMP1-TMP2)/(1-TMP2)];
    end
    Ratio_20_4_mean_ts = [Ratio_20_4_mean_ts;mean(tmp)]; 
    Ratio_20_4_std_ts = [Ratio_20_4_std_ts;std(tmp)];
end
%%
Ratio_21_4 = Ratio_21_2./Ratio_21_1;
Ratio_21_4_mean_ts = []; Ratio_21_4_std_ts = [];

for i = 0:5;
    tmp = [];
    for j = 1:3;
        if i ==0; jj=1;else;jj=j;end
        TMP1 = Ratio_21_1(find(Subset_21_1(:,5)==i & Subset_21_1(:,8)==jj));
        TMP2 = Ratio_21_2(find(Subset_21_2(:,5)==i & Subset_21_2(:,8)==jj));
        tmp = [tmp;(TMP1-TMP2)/(1-TMP2)];
    end
    Ratio_21_4_mean_ts = [Ratio_21_4_mean_ts;mean(tmp)]; 
    Ratio_21_4_std_ts = [Ratio_21_4_std_ts;std(tmp)];
end