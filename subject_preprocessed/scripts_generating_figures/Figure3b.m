load('S04_metagenomes/newstartend.mat');
sttedd=newsttedd;
load('S04_metagenomes/coveragematrix_meta.mat');
%%
Meta=[];
for i = 1:size(all_coverage_per_bp,1);
    i
    lowerquantile=quantile(all_coverage_per_bp(i,:),0.3);
    upperquantile=quantile(all_coverage_per_bp(i,:),0.7);
    metamean = mean(all_coverage_per_bp(i,find(all_coverage_per_bp(i,:)>lowerquantile & all_coverage_per_bp(i,:)<upperquantile)));
    
    m1 = mean(all_coverage_per_bp(i,sttedd(2,1):sttedd(2,2)));
    m2 = mean(all_coverage_per_bp(i,sttedd(3,1):sttedd(3,2)));
    Meta=[Meta;[m1/metamean m2/metamean]];
end
%%
load('S04_metagenomes/coveragematrix.mat');
%%
ISO=[];
for i = 1:size(all_coverage_per_bp,1);
    i
    isomean=mean(all_coverage_per_bp(i,:));
    m1 = mean(all_coverage_per_bp(i,sttedd(2,1):sttedd(2,2)));
    m2 = mean(all_coverage_per_bp(i,sttedd(3,1):sttedd(3,2)));
    ISO=[ISO;[m1/isomean m2/isomean]];
end


%% Plot anthor distribution of depth
Blue=[192 211 216]/255; Red=[212 160 156]/255; blk=[0.1 0.1 0.1];
Colors=[Blue;Red];
%% Plot two genetic regions that are putative mobile elements; also the two regions are very likely to be linked

figure(1);hold on;
subplot(1,2,1);hold on
MS=30;
orders=[ones(24,1)*1;ones(22,1)*2];
a=1:46;
reorder=a(randperm(length(a)));
jig=0.15
for j = 1:46;
    i=reorder(j)
    xax = 1+0.6*jig*randn(1);
    scatter(xax,log(max(ISO(i,1),0.01)),MS,Colors(orders(i),:),'filled','MarkerFaceAlpha',6/8);
end
% for i = 25:46;
%     xax = 1+jig*randn(1);
%     scatter(xax,log(max(ISO(i,1),0.01)),MS,Red,'filled','MarkerFaceAlpha',6/8);
% end
for i = 1;
    xax = 2+jig*randn(1);
    plot(xax,log(Meta(i,1)),'s','MarkerFaceColor',Blue,'MarkerEdgeColor',Blue,'MarkerSize',MS/3);
end
for i = 2;
    xax = 2+0.6*jig*randn(1);
    plot(xax,log(Meta(i,1)),'s','MarkerFaceColor',Red,'MarkerEdgeColor',Red,'MarkerSize',MS/3);
end

set(gca,'ytick',[log(0.01) log(0.1) log(1) log(10) log(100)]);
set(gca,'Yticklabel',{ '<0.01' '0.1'  '1' '10' '100'})
ylim([log(0.005) log(105)])
xlim([0.5 2.5])
set(gca,'xtick',[]);
set(gca,'TickDir','out')
set(gca,'FontSize',12,'LineWidth',1)
%
subplot(1,2,2);hold on
for j = 1:46;
    i=reorder(j)
    xax = 1+0.6*jig*randn(1);
    scatter(xax,log(max(ISO(i,2),0.01)),MS,Colors(orders(i),:),'filled','MarkerFaceAlpha',6/8);
end
for i = 1;
    xax = 2+jig*randn(1);
    plot(xax,log(Meta(i,2)),'s','MarkerFaceColor',Blue,'MarkerEdgeColor',Blue,'MarkerSize',MS/3);
end
for i = 2;
    xax = 2+jig*randn(1);
    plot(xax,log(Meta(i,2)),'s','MarkerFaceColor',Red,'MarkerEdgeColor',Red,'MarkerSize',MS/3);
end

set(gca,'ytick',[log(0.01) log(0.1) log(1) log(10) log(100)]);
set(gca,'Yticklabel',{ '<0.01' '0.1'  '1' '10' '100'})
ylim([log(0.005) log(105)])
xlim([0.5 2.5])
set(gca,'xtick',[]);
set(gca,'TickDir','out')
set(gca,'FontSize',12,'LineWidth',1)