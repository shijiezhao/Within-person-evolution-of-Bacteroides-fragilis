cd('S01_metagenomes')
%% Step 1: load data and scripts
load('coveragematrix.mat');
load('candidate_mutation_table.mat');
load('TimeSeries.mat');
load('optc.mat');   % Include the information of ancestor allele
load('Isolate_Info.mat');
load('change.mat');
load('IsolateMatrix.mat');
load('BF.mat')
load('TC.mat')
%% Initialization
counts=counts(:,change,:);
Quals=Quals(change,:);p=p(change);
MUT=[];TOT=[];N=-1;
bin=1;smt=0;mk=0;MaxD=0;

TOT7=[];MUT7=[];
for i = 1:82;
    [MIC,Days,MutC,TotC] = GetSNPtrajectory(i,counts,bin,N,Date,optc,smt);
    label = optc(i,4); % Where the mutation is
    if label == 2;
        DS=Days;CS=MIC;
        MUT=[MUT;MutC];
        TOT=[TOT;TotC];
    end
    if label == 3;
        label
        MUT=[MUT;TotC-MutC];
        TOT=[TOT;TotC];
    end
end
%%
QualityDates=[];
S1S2=[];
for i = 1:206;
    if sum(TOT(:,1)>20);
        QualityDates=[QualityDates;Days(i)];
        S1S2=[S1S2;double(sum(MUT(:,i)))/(sum(TOT(:,i))-sum(MUT(:,i)))];
    end
end

%%  The frequency of B.fragilis overtime
figure(1);
subplot(3,1,1);hold on;
lw=1.5
Blue=[192 211 216]/255; Red=[212 160 156]/255; blk=[0.2 0.2 0.2]*1;
MS=20;
plot([-20 600],[mean(BF) mean(BF)],'-','Color',blk*3,'Linewidth',lw);
plot(Days,BF,'-','Color',Blue,'Linewidth',1);
plot(Days,BF,'.','MarkerFaceColor',blk,'MarkerEdgeColor',blk,'MarkerSize',MS);
box off
set(gca,'ytick',[0 0.05 0.1]);set(gca,'Yticklabel',{'0' '5%' '10%'})
set(gca,'Xtick',[0 100 200 300 400 500 600]);set(gca,'Xticklabel',{'0' '100' '200' '300' '400' '500' '600'})
ylim([0 0.1]);xlim([-20 600]);set(gca,'TickDir','out');set(gca,'FontSize',12,'LineWidth',1)

%% The ratio of sl1/sl2 overtime
subplot(3,1,2);hold on;
MS=20;
TSM = 1./S1S2;
TSM=TSM(TSM<1000000);
plot([-20 600],[ mean(log(TSM)) mean(log(TSM))],'-','Color',blk*3,'Linewidth',lw);


plot(QualityDates,log(1./S1S2),'-','Color',Blue,'Linewidth',1);
plot(QualityDates,log(1./S1S2),'.','MarkerFaceColor',blk,'MarkerEdgeColor',blk,'MarkerSize',MS);
box off

set(gca,'ytick',[log(0.1) log(1) log(10) log(100)]);set(gca,'Yticklabel',{'0.1' '1' '10' '100'})
set(gca,'Xtick',[0 100 200 300 400 500 600]);set(gca,'Xticklabel',{'0' '100' '200' '300' '400' '500' '600'})

ylim([log(0.1) log(100)]);xlim([-20 600]);set(gca,'TickDir','out');set(gca,'FontSize',12,'LineWidth',1)

%% Plot D52
cd('../S03_metagenomes')
%% Step 1: load data and scripts
load('coveragematrix.mat');load('BfragAbundance.mat');load('TimeSeries.mat');

TotCounts=[];
for i = 1:74;
    i
    tmp = all_coverage_per_bp(i,:);
    mm = median(tmp);
    ms = sum(tmp(tmp<5*mm));    % regions that are larger than 5*median are probably from mobile element sequence, and are thus from other species
    TotCounts=[TotCounts;ms];
end
%%
subplot(3,1,3);hold on;
%
% TS1=Meta_counts(:,1)./Meta_counts(:,2);
TS2=(TotCounts/200)./Meta_counts(:,2);
plot(Date,TS2)

sep1=230;sep2=370;
exd=10;
plot([-5 150],[mean(TS2) mean(TS2)],'-','Color',blk*3,'Linewidth',lw);

plot(Date,TS2,'-','Color',Blue,'Linewidth',1);
plot(Date,TS2,'.','MarkerFaceColor',blk,'MarkerEdgeColor',blk,'MarkerSize',MS);

box off

set(gca,'ytick',[0 0.10 0.20]);set(gca,'Yticklabel',{'0'  '10%'  '20%'})

set(gca,'Xtick',[0 50 100 150]);set(gca,'Xticklabel',{'0' '50' '100' '150' })
xlim([-5 150]);ylim([0 0.2]);set(gca,'TickDir','out');set(gca,'FontSize',12,'LineWidth',1)