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
figure(1);hold on;
NormOtherLineages=zeros(6,11);
TMP1=[0;0;0;0;0;0]; % KEEP track of mutation spectrum for all other 11 lineages
for i = 1:11;
    tmp = [OtherLineages(1:3,i);OtherLineages(4:6,i)];
    TMP1 = TMP1+tmp;
    NormOtherLineages(:,i) = tmp/sum(tmp);
end
% Other lineage
c1=[1 1 1]*0.6;c0=[1 1 1]*0.2;
for i = 1:6;
    h=bar(i,mean(NormOtherLineages(i,:)),0.19);set(h,'FaceColor',c1);set(h,'LineWidth',0.05);set(h,'EdgeAlpha',0);
    h=plot([i,i],[mean(NormOtherLineages(i,:))-std(NormOtherLineages(i,:)),mean(NormOtherLineages(i,:))+std(NormOtherLineages(i,:))],'-','Color', c0, 'Linewidth',1)
end
% Normal lineage
c2=[234 201 157]/256;
tmp = [Normalsublineage(1:3);Normalsublineage(4:6)];
TMP2=tmp;
NormNorm=tmp/(sum(tmp));
for i = 1:6;
    h=bar(i +0.19,mean(NormNorm(i,:)),0.19);set(h,'FaceColor',c2);set(h,'LineWidth',0.05);set(h,'EdgeAlpha',0);
end
% Hyper lineage
c3=[139 112 165]/256;
tmp = [Hypersublineage(1:3);Hypersublineage(4:6)];
TMP3=tmp;
NormHyper=tmp/(sum(tmp));
for i = 1:6;
    xl = i;
    h=bar(xl+0.38,mean(NormHyper(i,:)),0.19);set(h,'FaceColor',c3);set(h,'LineWidth',0.05);set(h,'EdgeAlpha',0);
end

% Format the figure
xlim([0.5,7]);ylim([0 1]);set(gca,'YTick', [0 0.5 1]);set(gca,'XTick', {});set(gca,'XTickLabel',{});set(gca,'Fontsize',10);set( gca, 'TickDir', 'out' );set(gca,'linewidth',1)

% Calculate Chi-squre p-value
ChiSqPval(TMP1,TMP2,4) % compare normal branch and the combination of 11 other lineages
ChiSqPval(TMP1,TMP3,4) % compare hyper branch and the combination of 11 other lineages
