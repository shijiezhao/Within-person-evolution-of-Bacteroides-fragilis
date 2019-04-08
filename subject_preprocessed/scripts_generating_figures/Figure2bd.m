%%  Load files
load('summary_data/MEDi.mat') 
% In the data structure of MEDi:
% 1st column: lineage ID
% 2nd column: MED length
% 3rd column: present in % of isolates
% 4,5,6: gain,loss required to explain the phylogeny and inferred gain/loss
% based parsimony

%Initilization
SNPs=[];dMRCAT0=[];T2RS={}
for i = 1:12;
    if i < 10;
        load(['L0' num2str(i) '/SNP_dMRCA_related_intermediary.mat']);
    else
        load(['L' num2str(i) '/SNP_dMRCA_related_intermediary.mat']);
    end
    SNPs=[SNPs;nsnp];
    dMRCAT0=[dMRCAT0;dmrca_overtime(1)];
    T2RS{i}=Tip2Root;
end

%% Plot figure 1d, summarize dMRCA at initial sampling
lw=0.75; mycolor=[0.7 0.7 0.7];
figure(2);subplot(2,1,2);hold on;
for i = 1:12;
    if i~=8;
    h=bar(i,dMRCAT0(i),0.6);set(h,'FaceColor',mycolor);set(h,'LineWidth',lw);set(h,'EdgeAlpha',0.5);
    end
end
% Plot L08 TOTAL
% h=bar(8,dMRCAT0(8),0.5);set(h,'FaceColor',[1 1 1]*0.95);set(h,'LineWidth',lw);set(h,'EdgeAlpha',0.5);
% Plot L08, using only the normal sublineages
h=bar(8,9.9,0.5);set(h,'FaceColor',mycolor);set(h,'LineWidth',lw);set(h,'EdgeAlpha',0.5);

% format:
xlim([-0.25 12.75]);ylim([0 20]);set(gca,'Ytick',[0  2 4  6  8  10 14]);set(gca,'Yticklabel',{'0' '2'  '4' '6' '8' '10' '40'})
set(gca,'TickDir','out');set(gca,'Xtick',[0.25 12.75]);set(gca,'Xticklabel',{''  ''})

% 
C3=[1 1 1]*0.1;Colors={[51 102 153]/256; [164 164 164]/256;[175 65 27]/256;[25 25 25]/256}
jig1=0.08;
jig2=0.04
sz=30;
alp=.5;
alp1=.5;

%
% Plot figure 1d_NEW, summarize dMRCA at initial sampling, with boxplots + dots
lw=0.25; mycolor=[0.7 0.7 0.7];
figure(2);
subplot(2,1,1);hold on;subplot(2,1,2);hold on;
for i = 1:12;
    TIP2ROOTS=T2RS{i}{1};
    mth = mean(TIP2ROOTS(TIP2ROOTS<50));
    for j = 1:length(TIP2ROOTS);
        if TIP2ROOTS(j)<50;
            subplot(2,1,2);
            scatter(i+jig1*randn(1),TIP2ROOTS(j)+jig2*randn(1),sz,'MarkerEdgeColor',C3,'MarkerFaceColor',C3,'MarkerFaceAlpha',alp1,'MarkerEdgeAlpha',0)
        else
            subplot(2,1,1);
            scatter(i+jig1*randn(1),TIP2ROOTS(j)+jig2*randn(1),sz,'MarkerEdgeColor',Colors{3},'MarkerFaceColor',Colors{3},'MarkerFaceAlpha',alp1,'MarkerEdgeAlpha',0)
        end
    end
    subplot(2,1,2);
    plot([i-0.3,i+0.3],[mth,mth],'-','Color',Colors{3},'LineWidth',2)

end

% format:
subplot(2,1,2);
xlim([-0.25 12.75]);ylim([0 20])
set(gca,'TickDir','out');
subplot(2,1,1);
xlim([-0.25 12.75]);ylim([50 80]);
set(gca,'TickDir','out');


%2
half_wid = 0.25;
c1 = [0 0 0];
mean_length=1.5
for i = 1:12;
    TIP2ROOTS=T2RS{i}{1};Meta=TIP2ROOTS(TIP2ROOTS<50);
    ath = quantile(Meta,0.25);
    bth = quantile(Meta,0.5);
    cth = quantile(Meta,0.75);
    mth = mean(Meta)
    iqr = cth-ath;
    xl = i;
    plot([xl-half_wid,xl+half_wid],[ath,ath],'-','Color',c1,'LineWidth',1)
    plot([xl-half_wid,xl+half_wid],[bth,bth],'-','Color',c1,'LineWidth',1)
    plot([xl-half_wid,xl+half_wid],[cth,cth],'-','Color',c1,'LineWidth',1)
    plot([xl-half_wid,xl-half_wid],[ath,cth],'-','Color',c1,'LineWidth',1)
    plot([xl+half_wid,xl+half_wid],[ath,cth],'-','Color',c1,'LineWidth',1)
    plot([xl,xl],[ath,ath-1.5*iqr],'-','Color',c1,'LineWidth',1)
    plot([xl,xl],[cth,cth+1.5*iqr],'-','Color',c1,'LineWidth',1)
    plot([xl-half_wid/3,xl+half_wid/3],[ath-1.5*iqr,ath-1.5*iqr],'-','Color',c1,'LineWidth',1)
    plot([xl-half_wid/3,xl+half_wid/3],[cth+1.5*iqr,cth+1.5*iqr],'-','Color',c1,'LineWidth',1)
%     plot([xl-mean_length*half_wid,xl+mean_length*half_wid],[mth,mth],'-','Color',Colors{3},'LineWidth',2)
end

%% Plot figure 1b, summarize SNPs identified from each subject
figure(1);hold on;
subplot(2,1,1);hold on; % SNP
for i = 1:12;
    h=bar(i,SNPs(i),0.5); set(h,'FaceColor',mycolor);
    set(h,'LineWidth',lw); set(h,'EdgeAlpha',1);
end
% Plot L08
h=bar(8,120,0.5);set(h,'FaceColor',mycolor);set(h,'LineWidth',lw);set(h,'EdgeAlpha',1);
    
xlim([-0.25 12.75]);ylim([0 120]);set(gca,'Ytick',[0 20 40 60 80 100 120])
set(gca,'Yticklabel',{'0' '20'  '40' '60' '80' '100' num2str(max(SNPs))});
set(gca,'TickDir','out')
set(gca,'Xtick',[0.25 12.75])
set(gca,'Xticklabel',{''  ''})

%% Plot figure 1b, summarize MEDs identified from each subject
subplot(2,1,2);hold on;
c1=[100 142 179]/255;c2=[0.9 0.9 0.9];c3=[188 56 50]/255;
MEDs=zeros(12,3);
for i = 1:12;
    MEDs(i,1)=sum(MEDi(:,1)==i & MEDi(:,6)==1);     % gain
    MEDs(i,2)=sum(MEDi(:,1)==i & MEDi(:,6)==-1);    % loss
    MEDs(i,3)=sum(MEDi(:,1)==i & MEDi(:,6)==0);     % unclear
end
for i = 1:12;
    % Gain
    h=bar(i,sum(MEDs(i,:)),0.5);set(h,'FaceColor',c2);
    set(h,'LineWidth',lw); set(h,'EdgeAlpha',1);
    heit = sum(MEDs(i,:));
    for yv = 0:0.25:heit+0.5;
        x1 = max(0,yv-heit)+i-0.25; x2 = min(yv-0.5,0)+i+0.25;
        y1 = min(yv,heit); y2 = max(yv-0.5,0);
        plot([x1 x2],[y1 y2],'-k');

    end
    % Unclear
    h=bar(i,sum(MEDs(i,2:3)),0.5);set(h,'FaceColor',c2);set(h,'LineWidth',lw);set(h,'EdgeAlpha',1);
    % Loss
    h=bar(i,sum(MEDs(i,2)),0.5);set(h,'FaceColor',mycolor);set(h,'LineWidth',lw);set(h,'EdgeAlpha',1);
end
%% Format
xlim([-0.25 12.75]);ylim([0 9])
set(gca,'Ytick',[0 3 6 9])
set(gca,'Yticklabel',{'0' '3'  '6' '9'})
set(gca,'Xtick',[])
set(gca,'TickDir','out')
set(gca,'Xtick',[0.25 12.75])
set(gca,'Xticklabel',{''  ''})



