%% Load data
load('summary_data/GenomeInformationForMembraneEnrichmentAnalysis.mat')
% Data: 12 lineage genome, and the cellular localization prediction of the
% gene, 1: membrane related; 0: cytoplasmic; 2: extracellular
Simulation_Results=[];
totalmut=46; % Total number of mutations in the 16 genes
for Rpts = 1:1000; % 1000 simulation
    AllRandomGene=[];
    MembraneRelated=0;
    for gene = 1:16;
        genome = MutationTable{gene,1}(randi(length(MutationTable{gene,1}))); % randomly select a genome from the mutated gene
        % randomly pick a gene from that genome
        tmp = randi(length(Data{genome}));
        selectedGene=Data{genome}(tmp,1);
        localization=Data{genome}(tmp,2);
        if localization == 1;
            MembraneRelated=MembraneRelated+MutationTable{gene,2}; % Count the number of how many times this original gene was mutated
        end
    end
    Simulation_Results=[Simulation_Results;MembraneRelated/totalmut]; 
end

%%

Red1=[244 113 126]/256; Red2=[212 160 156]/256;
Blue1=[64 134 153]/256; Blue2=[193 219 216]/256; Blue3=[100 142 180]/256;
Grey1=[64 64 64]/256; Grey2=[200 200 200]/256;

c1 = Grey1; c2=Grey2;c3=Blue2;c4=Red2;
mycolors=([c1;c2;c3;c4]);



figure(1);hold on;
hold on;
pts=0:0.025:1;
[bincounts] = histc(Simulation_Results,pts);
h=bar(pts,bincounts,'histc')
h.FaceColor = c2;

% h = bar(39/49,300,0.01)
h = bar(36/46,300,0.01)

set(h,'FaceColor',Red2);
set(h,'EdgeAlpha',0);

xlim([0 1])
ylim([0 150])
set(gca,'Ytick',[0 50 100 150 200 250])
set(gca,'Xtick',[0 0.2 0.4 0.6 0.8 1.0])
set(gca,'Yticklabel',{})
set(gca,'Xticklabel',{})

set(gca,'TickDir','out')
set(gca,'FontSize',10)

saveas(gcf,['TotWithin.fig']);
set(gca,'linewidth',1)
