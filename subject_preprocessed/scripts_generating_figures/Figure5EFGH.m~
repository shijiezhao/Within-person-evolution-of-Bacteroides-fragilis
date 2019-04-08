load('competition_exp/Phage_competition.mat');
tmP = [0;9;18;27;36];
%% Perc = s1b%
Equals_avg=[];
Equals_err=[];
% figure(3);hold on;
for cmb = 1:18;
Subset_1 = Quantify2(Quantify2(:,6)==cmb,:);
Ratio_1 = (Subset_1(:,1)+Subset_1(:,2))./(Subset_1(:,1)+Subset_1(:,2)+Subset_1(:,3)+Subset_1(:,4));
Ratio_1_mean_ts = []; Ratio_1_std_ts = [];
for i = 0:4;
    Ratio_1_mean_ts = [Ratio_1_mean_ts; mean(Ratio_1(find(Subset_1(:,5)==i)))]; 
    Ratio_1_std_ts = [Ratio_1_std_ts; var(Ratio_1(find(Subset_1(:,5)==i)))];
end
Equals_avg=[Equals_avg; 1-Ratio_1_mean_ts'];
Equals_err=[Equals_err; Ratio_1_std_ts'];
% e=errorbar(tmP,Equals_avg(cmb,:),Equals_err(cmb,:),'LineWidth',1,'CapSize',10);
% pause
end
%%
figure(3);
subplot(1,3,3);hold on;
for i = 1:8;
    e=errorbar(tmP,Equals_avg(i,:),Equals_err(i,:),'LineWidth',1,'CapSize',10);
    e.Color = C1;
    e=errorbar(tmP,1-Equals_avg(i,:),Equals_err(i,:),'LineWidth',1,'CapSize',10);
    e.Color = C2;
end
xlim([0 36]);
ylim([0 1]);
plot([18 18],[0 1],'k--','LineWidth',1)
box on;
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'LineWidth',1)

%%
figure(3);
subplot(1,3,1);hold on;
for i = 9:16;
    e=errorbar(tmP,Equals_avg(i,:),Equals_err(i,:),'LineWidth',1,'CapSize',10);
    e.Color = C1;
    e=errorbar(tmP,1-Equals_avg(i,:),Equals_err(i,:),'LineWidth',1,'CapSize',10);
    e.Color = C2;
end
xlim([0 36]);
ylim([0 1]);
plot([18 18],[0 1],'k--','LineWidth',1)
box on; 
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'LineWidth',1)


%%
figure(3);
subplot(1,3,2);hold on;
for i = 17:18;
    e=errorbar(tmP,1-Equals_avg(i,:),Equals_err(i,:),'LineWidth',1,'CapSize',10);
    e.Color = C1;
    e=errorbar(tmP,Equals_avg(i,:),Equals_err(i,:),'LineWidth',1,'CapSize',10);
    e.Color = C2;
end
xlim([0 36]);
ylim([0 1]);
plot([18 18],[0 1],'k--','LineWidth',1)

box on; 
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'LineWidth',1)
