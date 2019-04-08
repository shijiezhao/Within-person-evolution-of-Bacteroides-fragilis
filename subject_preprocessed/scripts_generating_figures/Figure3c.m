load('S01_metagenomes/coveragematrix_meta.mat');
%% Metagenome coverage
load('S01_metagenomes/sttedd.mat')
Meta=[];
for i = 1:size(all_coverage_per_bp,1);
    i
    lowerquantile=quantile(all_coverage_per_bp(i,:),0.3);
    upperquantile=quantile(all_coverage_per_bp(i,:),0.7);
    metamean = mean(all_coverage_per_bp(i,find(all_coverage_per_bp(i,:)>lowerquantile & all_coverage_per_bp(i,:)<upperquantile)));
    
    m1 = mean(all_coverage_per_bp(i,sttedd(1,1):sttedd(1,2)));
    m2 = mean(all_coverage_per_bp(i,sttedd(2,1):sttedd(2,2)));
    Meta=[Meta;[m1/metamean m2/metamean]];
end
%% Isolates coverage
load('L01/coveragematrix.mat');
%%
ISO=[];
for i = 1:size(all_coverage_per_bp,1);
    i
    isomean=mean(all_coverage_per_bp(i,:));
    m1 = mean(all_coverage_per_bp(i,sttedd(1,1):sttedd(1,2)));
    m2 = mean(all_coverage_per_bp(i,sttedd(2,1):sttedd(2,2)));
    ISO=[ISO;[m1/isomean m2/isomean]];
end

%% Plot anthor distribution of depth
Purple1=[138 112 165]/255; Yellow=[233 202 163]/255; blk=[0.1 0.1 0.1];
Grey=[1 1 1]*0.7;

figure(1);hold on;
MS=40;

jig=0.2
for i = 1:187;
    xax = 1+0.6*jig*randn(1);
    scatter(xax,log(ISO(i,1)/2+ISO(i,2)/2),MS,Grey,'filled','MarkerFaceAlpha',3/8);
end
for i = 1:206;
    xax = 2+0.6*jig*randn(1);
    scatter(xax,log(Meta(i,1)/2+Meta(i,2)/2),MS,Grey,'filled','MarkerFaceAlpha',3/8);
end
for i = 1:415;
    xax = 3+0.6*jig*randn(1);
    scatter(xax,log(0.01),MS,Grey,'filled','MarkerFaceAlpha',3/8);
end

set(gca,'ytick',[log(0.01) log(0.1) log(1) log(10) log(100) log(300)]);
set(gca,'Yticklabel',{ '<0.01' '0.1'  '1' '10' '100' '300'})
ylim([log(0.005) log(300)])
xlim([0.5 3.5])
set(gca,'xtick',[]);
set(gca,'TickDir','out')
set(gca,'FontSize',12,'LineWidth',1)
