%% Read in files to plot (a-c)
L01ts=[];L01mdSNPs=[];ts=[];mdSNPs=[];
for i = 1; % Read L01 average
    load(['L0' num2str(i) '/SNP_dMRCA_related_intermediary.mat']);
    for ltp = 2:length(t);
        L01ts = [L01ts t(ltp)/365];
        L01mdSNPs=[L01mdSNPs mean(dSNP{ltp})];
    end
end
for i = 2:7;
        load(['L0' num2str(i) '/SNP_dMRCA_related_intermediary.mat']);
    for ltp = 2:length(t);
        ts = [ts t(ltp)/365];
        mdSNPs=[mdSNPs mean(dSNP{ltp})];
    end
end
%% Keep record the linear fit for panel a-h
Rsqures=[]; LinearFits=[]; 
C3=[1 1 1]*0.6;C1=[1 1 1]*0.35;
%% For supplementary figure 5; Panel (a). calculate molecular clock using average dSNPs from every time
subplot(4,4,1);hold on; 
% point
x=[L01ts ts];y=[L01mdSNPs mdSNPs];
[B,S] = polyfit(x,y, 1); p = polyfit(x, y, 1); yfit = polyval(p,x); yresid = y - yfit; SSresid = sum(yresid.^2); SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal;
Rsqures=[Rsqures;rsq];LinearFits=[LinearFits;p];
% Plot the fitting curve first
plot([0 2],[0*B(1)+B(2) 2*B(1)+B(2)],'-', 'Color', C3, 'Linewidth',1);
% for L01:
h=plot(L01ts,L01mdSNPs,'o','Color',C1,'MarkerSize',5); set(h,'MarkerEdgeColor',C1,'MarkerFaceColor',C1);
% for non-L01:
h=plot(ts,mdSNPs,'o','Color',C2,'MarkerSize',5); set(h,'MarkerEdgeColor',C2,'MarkerFaceColor',C2);
% Format the figure
xlim([0,2]); ylim([0 2]); set(gca,'YTick', [0 1 2]);set(gca,'Fontsize',10);set( gca, 'TickDir', 'out' );set(gca,'linewidth',1);

%% Panel (b). calculate molecular clock using average dSNPs using time points from L01 only
subplot(4,4,2);hold on; % L01-only; average
% Fitting!
x=[L01ts];y=[L01mdSNPs];
[B,S] = polyfit(x,y, 1); p = polyfit(x, y, 1); yfit = polyval(p,x); yresid = y - yfit; SSresid = sum(yresid.^2); SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal;
Rsqures=[Rsqures;rsq];LinearFits=[LinearFits;p];
% Plot the fitting curve first
plot([0 2],[0*B(1)+B(2) 2*B(1)+B(2)],'-', 'Color', C3, 'Linewidth',1);
% for L01:
h=plot(L01ts,L01mdSNPs,'o','Color',C1,'MarkerSize',5); set(h,'MarkerEdgeColor',C1,'MarkerFaceColor',C1);
% Format the figure
xlim([0,2]); ylim([0 2]); set(gca,'YTick', [0 1 2]);set(gca,'Fontsize',10);set( gca, 'TickDir', 'out' );set(gca,'linewidth',1);

%% Panel (c). calculate molecular clock using average dSNPs using time points from non-L01 lineages
subplot(4,4,3);hold on; % L01-only; average
% Fitting!
x=[ts];y=[mdSNPs];
[B,S] = polyfit(x,y, 1); p = polyfit(x, y, 1); yfit = polyval(p,x); yresid = y - yfit; SSresid = sum(yresid.^2); SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal;
Rsqures=[Rsqures;rsq];LinearFits=[LinearFits;p];
% Plot the fitting curve first
plot([0 2],[0*B(1)+B(2) 2*B(1)+B(2)],'-', 'Color', C3, 'Linewidth',1);
% for non-L01:
h=plot(ts,mdSNPs,'o','Color',C2,'MarkerSize',5); set(h,'MarkerEdgeColor',C2,'MarkerFaceColor',C2);
% Format the figure
xlim([0,2]); ylim([0 2]); set(gca,'YTick', [0 1 2]);set(gca,'Fontsize',10);set( gca, 'TickDir', 'out' );set(gca,'linewidth',1);


%% Panel (d-f). calculate molecular clock using all isolate's dSNPs without averages
ytop=7; % define a y limit
ts=[];idSNPs=[];L01ts=[];L01idSNPs=[]; 
for i = 1:12;
    if i<10;
    	load(['L0' num2str(i) '/SNP_dMRCA_related_intermediary.mat']);
    else;
        load(['L' num2str(i) '/SNP_dMRCA_related_intermediary.mat']);
    end
    for ltp = 2:length(t);
        for j = 1:length(dSNP{ltp});
            if i==1; % for L01
                L01ts = [L01ts t(ltp)/365];
                L01idSNPs=[L01idSNPs dSNP{ltp}(j)];
            else;    % for non-L01
                ts = [ts t(ltp)/365];
                idSNPs=[idSNPs dSNP{ltp}(j)];
            end
        end
    end
end
jig=0.15; % Many isolates from a same time point have same dSNP values, give them jigs to visualize

%% Panel (d). calculate molecular clock using all isolate's dSNPs from all time points from all donors
subplot(4,4,4);hold on; % All; individual
% Fitting!
x=[L01ts ts];y=[L01idSNPs idSNPs];
[B,S] = polyfit(x,y, 1); p = polyfit(x, y, 1); yfit = polyval(p,x); yresid = y - yfit; SSresid = sum(yresid.^2); SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal;
Rsqures=[Rsqures;rsq];LinearFits=[LinearFits;p];
% Plot the fitting curve first
plot([0 2],[0*B(1)+B(2) 2*B(1)+B(2)],'-', 'Color', C3, 'Linewidth',1);
% for L01:
for i = 1:length(L01ts);
    h=plot(L01ts(i)+jig*rand(1)/3,L01idSNPs(i)+jig*rand(1),'o','Color',C1,'MarkerSize',2);
    set(h,'MarkerEdgeColor',C1,'MarkerFaceColor',C1);
end
% for non-L01:
for i = 1:length(ts);
    h=plot(ts(i)+jig*rand(1)/3,idSNPs(i)+jig*rand(1),'o','Color',C2,'MarkerSize',2);
    set(h,'MarkerEdgeColor',C2,'MarkerFaceColor',C2);
end
% Format the figure
xlim([0,2]);ylim([0 ytop]);set(gca,'Fontsize',10);set( gca, 'TickDir', 'out' );set(gca,'linewidth',1)

%% Panel (e). calculate molecular clock using all isolate's dSNPs from L01 only
subplot(4,4,5);hold on; % L01; individual
% Fitting!
x=[L01ts];y=[L01idSNPs];
[B,S] = polyfit(x,y, 1); p = polyfit(x, y, 1); yfit = polyval(p,x); yresid = y - yfit; SSresid = sum(yresid.^2); SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal;
Rsqures=[Rsqures;rsq];LinearFits=[LinearFits;p];
% Plot the fitting curve first
plot([0 2],[0*B(1)+B(2) 2*B(1)+B(2)],'-', 'Color', C3, 'Linewidth',1);
% for L01:
for i = 1:length(L01ts);
    h=plot(L01ts(i)+jig*rand(1)/3,L01idSNPs(i)+jig*rand(1),'o','Color',C1,'MarkerSize',2);
    set(h,'MarkerEdgeColor',C1,'MarkerFaceColor',C1);
end
% Format the figure
xlim([0,2]);ylim([0 ytop]);set(gca,'Fontsize',10);set( gca, 'TickDir', 'out' );set(gca,'linewidth',1)


%% Panel (f). calculate molecular clock using all isolate's dSNPs from non-L01 lineages
subplot(4,4,6);hold on; % All; individual
% Fitting!
x=[ts];y=[idSNPs];
[B,S] = polyfit(x,y, 1); p = polyfit(x, y, 1); yfit = polyval(p,x); yresid = y - yfit; SSresid = sum(yresid.^2); SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal;
Rsqures=[Rsqures;rsq];LinearFits=[LinearFits;p];
% Plot the fitting curve first
plot([0 2],[0*B(1)+B(2) 2*B(1)+B(2)],'-', 'Color', C3, 'Linewidth',1);
% for non-L01:
for i = 1:length(ts);
    h=plot(ts(i)+jig*rand(1)/3,idSNPs(i)+jig*rand(1),'o','Color',C2,'MarkerSize',2);
    set(h,'MarkerEdgeColor',C2,'MarkerFaceColor',C2);
end
% Format the figure
xlim([0,2]);ylim([0 ytop]);set(gca,'Fontsize',10);set( gca, 'TickDir', 'out' );set(gca,'linewidth',1)

%% Panel (g), using tip to root over time to calculate molecular clock, using data from L01
ts=[];mT2Rs=[]; % Initiate
load(['L01/SNP_dMRCA_related_intermediary.mat']);
for ltp = 1:length(t);
    ts = [ts t(ltp)/365]; mT2Rs=[mT2Rs t2r_overtime(ltp)];
end
subplot(4,4,7);hold on; % All; individual
% Fitting!
x=[ts];y=[mT2Rs];
[B,S] = polyfit(x,y, 1); p = polyfit(x, y, 1); yfit = polyval(p,x); yresid = y - yfit; SSresid = sum(yresid.^2); SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal; Rsqures=[Rsqures;rsq];LinearFits=[LinearFits;p];
% Plot the fitting curve first
plot([0 2],[0*B(1)+B(2) 2*B(1)+B(2)],'-', 'Color', C3, 'Linewidth',1);
% for L01:
h=plot(ts,mT2Rs,'o','Color',C1,'MarkerSize',5);set(h,'MarkerEdgeColor',C1,'MarkerFaceColor',C1);
% Format the figure
xlim([0,2]);ylim([7 11]);set(gca,'YTick', [7 8 9 10 11]);set(gca,'XTick', [0 1 2]);set(gca,'Fontsize',10);set( gca, 'TickDir', 'out' );set(gca,'linewidth',1);

%% Panel (h), using tip to root over time to calculate molecular clock, using data from L05
ts=[];mT2Rs=[];load(['L05/SNP_dMRCA_related_intermediary.mat']);
for ltp = 1:length(t);
    ts = [ts t(ltp)/365]; mT2Rs=[mT2Rs t2r_overtime(ltp)];
end
subplot(4,4,8);hold on; % All; individual
% Fitting!
x=[ts(1:3)];y=[mT2Rs(1:3)];
[B,S] = polyfit(x,y, 1); p = polyfit(x, y, 1); yfit = polyval(p,x); yresid = y - yfit; SSresid = sum(yresid.^2); SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal; Rsqures=[Rsqures;rsq];LinearFits=[LinearFits;p];
% Plot the fitting curve first
plot([0 2],[0*B(1)+B(2) 2*B(1)+B(2)],'-', 'Color', C3, 'Linewidth',1);
% for L05:
h=plot(x,y,'o','Color',C1,'MarkerSize',5);set(h,'MarkerEdgeColor',C1,'MarkerFaceColor',C1);
% Format the figure
xlim([0 1]);ylim([4 7]);set(gca,'YTick', [4 5 6 7 8]);set(gca,'XTick', [0 .5 1]);set(gca,'Fontsize',10);set( gca, 'TickDir', 'out' );set(gca,'linewidth',1)


%% Panel (i-o), dMRCA over time for L01-L07
for i = 1:7;
    load(['L0' num2str(i) '/SNP_dMRCA_related_intermediary.mat']);
    subplot(4,4,i+8);hold on;set(gca,'Fontsize',10);set(gca,'linewidth',1);set(gca,'TickDir','out');
    ty = t/365; % TIME IN YEAR
    for j = 1:size(dmrca_overtime,1);
        sem = dmrca_overtime(j,2)/sqrt(length(dSNP{j})); % Standard error of the mean
        h=plot([ty(j),ty(j)],[dmrca_overtime(j,1)-sem,dmrca_overtime(j,1)+sem]','-', 'Color', [.75 .75 .75], 'Linewidth',1);
    end
    h=plot(ty,dmrca_overtime(:,1),'-', 'Color', [.75 .75 .75], 'Linewidth',2);
    h=plot(ty,dmrca_overtime(:,1),'o','Color',[0.5 0.5 0.5],'MarkerSize',5);
    set(h,'MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor',[0.5 0.5 0.5]);
    xlim([-max(ty)*0.2 max(ty)*1.2]);
%     ylim([min(dmrca_overtime(:,1))-(max(dmrca_overtime(:,1))-min(dmrca_overtime(:,1)))*0.8 max(dmrca_overtime(:,1))+(max(dmrca_overtime(:,1))-min(dmrca_overtime(:,1)))*0.8]);
    ylim([min(dmrca_overtime(:,1)-dmrca_overtime(:,2)/4)-0.15,max(dmrca_overtime(:,1)+dmrca_overtime(:,2)/4)+0.15])
    box off;  
end

%% Panel (p): plot the mutational spectrum when excluding the GC-TA hypermutation type
% Create a dictionary
c=containers.Map; 
c('AT')=1;c('TA')=1;c('AC')=2;c('TG')=2;c('GC')=3;c('CG')=3; 
c('GT')=4;c('CA')=4;c('AG')=5;c('TC')=5;c('GA')=6;c('CT')=6;
OtherLineages=zeros(6,11); % mutation spectrum of non-L08 lineages
for i = 1:12;
    if i ~=8; % excluding L08
        if i <10;
            load(['L0' num2str(i) '/S0' num2str(i) 'prokka.mat']);
        else
            load(['L' num2str(i) '/S' num2str(i) 'prokka.mat']);
        end
        for j = 1:length(annotation_full);
            nt1 = annotation_full(j).anc; nts = annotation_full(j).nts;
            if nt1==nts(1);
                nt2=nts(2);
            else;
                nt2=nts(1);
            end
            type = c([nt1 nt2]); % determine mutation type
            if i < 8;
                OtherLineages(type,i)=OtherLineages(type,i)+1;
            else
                OtherLineages(type,i-1)=OtherLineages(type,i-1)+1;
            end
        end
    end
end
% Normal sublineage in L08
Normalsublineage=[0;0;0;0;0;0];
load(['L08/S08_norm.mat']);
for j = 1:length(annotation_full);
    nt1 = annotation_full(j).anc;
    nts = annotation_full(j).nts;
    if nt1==nts(1);
        nt2=nts(2);
    else;
        nt2=nts(1);
    end
    type = c([nt1 nt2]);
    Normalsublineage(type)=Normalsublineage(type)+1;
end
% Normal sublineage in L08
Hypersublineage=[0;0;0;0;0;0];
load(['L08/S08_hyper.mat']);
for j = 1:length(annotation_full);
    nt1 = annotation_full(j).anc;
    nts = annotation_full(j).nts;
    if nt1==nts(1);
        nt2=nts(2);
    else;
        nt2=nts(1);
    end
    type = c([nt1 nt2]);
    Hypersublineage(type)=Hypersublineage(type)+1;
end


%% Panel (p) continued: plot the mutational spectrum when excluding the GC-TA hypermutation type
subplot(4,4,16);hold on;
NormOtherLineages=zeros(5,11);
TMP1=[0;0;0;0;0]; % KEEP track of mutation spectrum for all other 11 lineages
for i = 1:11;
    tmp = [OtherLineages(1:3,i);OtherLineages(5:6,i)];
    TMP1 = TMP1+tmp;
    NormOtherLineages(:,i) = tmp/sum(tmp);
end
% Other lineage
c1=[1 1 1]*0.6;c0=[1 1 1]*0.2;
for i = 1:5;
    h=bar(i+(i>3),mean(NormOtherLineages(i,:)),0.19);set(h,'FaceColor',c1);set(h,'LineWidth',0.05);set(h,'EdgeAlpha',0);
    h=plot([i+(i>3),i+(i>3)],[mean(NormOtherLineages(i,:))-std(NormOtherLineages(i,:)),mean(NormOtherLineages(i,:))+std(NormOtherLineages(i,:))],'-','Color', c0, 'Linewidth',1)
end
% Normal lineage
c2=[234 201 157]/256;
tmp = [Normalsublineage(1:3);Normalsublineage(5:6)];
TMP2=tmp;
NormNorm=tmp/(sum(tmp));
for i = 1:5;
    h=bar(i + (i>3)+0.19,mean(NormNorm(i,:)),0.19);set(h,'FaceColor',c2);set(h,'LineWidth',0.05);set(h,'EdgeAlpha',0);
end
% Hyper lineage
c3=[139 112 165]/256;
tmp = [Hypersublineage(1:3);Hypersublineage(5:6)];
TMP3=tmp;
NormHyper=tmp/(sum(tmp));
for i = 1:5;
    xl = i + (i>3);
    h=bar(xl+0.38,mean(NormHyper(i,:)),0.19);set(h,'FaceColor',c3);set(h,'LineWidth',0.05);set(h,'EdgeAlpha',0);
end

% Format the figure
xlim([0.5,7]);ylim([0 1]);set(gca,'YTick', [0 0.5 1]);set(gca,'XTick', {});set(gca,'XTickLabel',{});set(gca,'Fontsize',10);set( gca, 'TickDir', 'out' );set(gca,'linewidth',1)

% Calculate Chi-squre p-value
ChiSqPval(TMP1,TMP2,4) % compare normal branch and the combination of 11 other lineages
ChiSqPval(TMP1,TMP3,4) % compare hyper branch and the combination of 11 other lineages





