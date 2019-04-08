dMRCAs={};antnts = ancnti(goodpos);
TempGoodsamples=goodsamples;
repeats=100;
figure(11);hold on;
set(gca,'LineWidth',1);
set( gca, 'TickDir', 'out' );
xlabel('Number of isolates')
ylabel('dMRCA')
SNPs={};
c1 = [231 146 110]/256;
c2 = [121 196 221]/256;
%% for each number from 1 to number of samples
for subsampling = 1:length(SampleNames(~in_outgroup));
    dMRCAs{subsampling} = []; SNPs{subsampling} = [];
    subsampling; % Number
    for r = 1:repeats; %numel(samplestoplot);
        dmrca = 0;
        antnts = ancnti(goodpos);
        Subsets = 1:length(SampleNames(~in_outgroup));
        ToRemoveSample=randsample(Subsets,length(Subsets)-subsampling);
        Subsets(ToRemoveSample)=0;
        Rsub = [];
        for tt = 1:length(AllSampNames);
            for abc = 1:length(Subsets);
                if Subsets(abc)>0;
                    if strcmp(SampleNames{Subsets(abc)},AllSampNames{tt});
                        Rsub = [Rsub;tt];
                    end
                end
            end
        end
        % Rsub: the set of isolates in the subset
        % length(Rub) = subsampling
        Indeces=[];
        for s = 1:subsampling;
            sample = Rsub(s);
            [val indx] = max(Allcounts(1:4,goodpos,sample) + Allcounts(5:8,goodpos,sample));
            Indeces = [Indeces;indx];
        end
        snpn=0;dmrca=0;
        for l = 1:length(antnts);
            if (min(Indeces(:,l)) ~= max(Indeces(:,l)));
                if antnts(l)==0;
                    antnts(l)=mode(Indeces(:,l));
                end
                dmrca = dmrca + sum(antnts(l) ~= Indeces(:,l))/subsampling;
                snpn = snpn+1;
            end
        end
        dMRCAs{subsampling}=[dMRCAs{subsampling};dmrca];
        SNPs{subsampling} = [SNPs{subsampling};snpn];
    end
    scatter(subsampling,mean(dMRCAs{subsampling}),50,'filled', ...
       'MarkerFaceAlpha',1,'MarkerFaceColor',c1);
    plot([subsampling subsampling],[mean(dMRCAs{subsampling})-std(dMRCAs{subsampling}) mean(dMRCAs{subsampling})+std(dMRCAs{subsampling})],'-','Color',[0.6,0.6,0.6],'LineWidth',1);
end
saveas(gcf,[Donor '_dMRCA.png']);
figure(12);hold on;
set(gcf, 'PaperSize', [4 2]);
set(gca,'LineWidth',1);
set( gca, 'TickDir', 'out' );
xlabel('Number of isolates')
ylabel('Number of SNPs')
for subsampling = 1:length(goodsamples(~in_outgroup));
    scatter(subsampling,mean(SNPs{subsampling}),50,'filled', ...
       'MarkerFaceAlpha',1,'MarkerFaceColor',c2);
    plot([subsampling subsampling],[mean(SNPs{subsampling})-std(SNPs{subsampling}) mean(SNPs{subsampling})+std(SNPs{subsampling})],'-','Color',[0.6,0.6,0.6],'LineWidth',1);
end
saveas(gcf,[Donor '_SNPs.png']);
save([Donor '_collector_curves.mat'],'dMRCAs','SNPs');
mkdir('Collector')
eval(['! mv ' Donor '* Collector'])