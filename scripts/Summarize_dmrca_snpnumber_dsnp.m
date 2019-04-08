
%% New Step 2: plot dMRCA vs #iso and SNP vs #iso for each time point
fignumdmrca = 123;
fignumsnp = 1024;
figure(fignumdmrca);hold on;
set(gca,'LineWidth',1);
set( gca, 'TickDir', 'out' );
xlabel('Number of isolates')
ylabel('dMRCA')
figure(fignumsnp);hold on;
set(gca,'LineWidth',1);
set( gca, 'TickDir', 'out' );
xlabel('Number of isolates')
ylabel('SNP')
c1 = [193 220 216]/256;
c2 = [212 160 156]/256;
c3 = [50 92 127]/256;
c4 = [234 201 157]/256;
c5 = [139 112 165]/256;
c6 = [64 64 64]/256;
c=[c1;c2;c3;c4;c5;c6];
dmrca_overtime=[];
stddmrca_overtime=[];
jig = 0.05
for ltp = 1:snum;
    dMRCAs={};antnts = ancnti(goodpos);
    stddMRCAs={};
    TempGoodsamples=TPS{ltp}';
    repeats=600;
    
    SNPs={};
    
    for subsampling = 1:sum(TempGoodsamples);
        dMRCAs{subsampling} = []; SNPs{subsampling} = [];stddMRCAs{subsampling}=[];
        for r = 1:repeats; %numel(samplestoplot);
            dmrca = 0;
            totdmrca=[];
            antnts = ancnti(goodpos);
            Subsets = find(TempGoodsamples>0);
            SelectedSample=randsample(Subsets,subsampling);
            SelMajorNT = majorNT(goodpos,SelectedSample);
            snpn=0;dmrca=0;
            for l = 1:length(antnts);
                if (min(SelMajorNT(l,:)) ~= max(SelMajorNT(l,:))); % Two alleles
%                     if antnts(l)==0;
%                         antnts(l)=mode(Indeces(:,l));
%                     end
                    dmrca = dmrca + sum(antnts(l) ~= SelMajorNT(l,:))/subsampling;
                    snpn = snpn+1;
                end
            end
            dMRCAs{subsampling}=[dMRCAs{subsampling};dmrca];
            SNPs{subsampling} = [SNPs{subsampling};snpn];
        end
        figure(fignumdmrca)
        scatter(subsampling+jig*ltp,mean(dMRCAs{subsampling}),50,'filled','MarkerFaceAlpha',1,'MarkerFaceColor',c(ltp,:));
        a = plot([subsampling+jig*ltp subsampling+jig*ltp],[mean(dMRCAs{subsampling})-std(dMRCAs{subsampling}) mean(dMRCAs{subsampling})+std(dMRCAs{subsampling})],'-','Color',c(ltp,:),'LineWidth',1);
        a.Color(4) = 0.2;
        figure(fignumsnp)
        scatter(subsampling+jig*ltp,mean(SNPs{subsampling}),50,'filled','MarkerFaceAlpha',1,'MarkerFaceColor',c(ltp,:));
        a = plot([subsampling+jig*ltp subsampling+jig*ltp],[mean(SNPs{subsampling})-std(SNPs{subsampling}) mean(SNPs{subsampling})+std(SNPs{subsampling})],'-','Color',c(ltp,:),'LineWidth',1);
        a.Color(4) = 0.2;
    end
    for s = 1:subsampling;
        
        totdmrca = [totdmrca;sum(antnts ~= SelMajorNT(:,s))];
    end
    dmrca_overtime = [dmrca_overtime; mean(totdmrca)];
    stddmrca_overtime = [stddmrca_overtime;std(totdmrca)];
end


if snum>1;
%% New step 3: plot the dMRCA over time
figure(545);hold on;
set(gca,'LineWidth',1);
set( gca, 'TickDir', 'out' );
ylabel('dMRCA')
for ltp = 1:snum;
    h=bar(ltp,dmrca_overtime(ltp),0.5);
    plot([ltp ltp],[dmrca_overtime(ltp)-stddmrca_overtime(ltp), dmrca_overtime(ltp)+stddmrca_overtime(ltp)],'-','Color',[0.3 0.3 0.3],'LineWidth',2);
    set(h,'FaceColor',c(ltp,:));
    set(h,'LineWidth',0.05);
    set(h,'EdgeAlpha',0);
end
dif=0.5
xlim([1-dif snum+dif])
ylim([0 3])
set(gca,'Ytick',[0 2.5 5 7.5 10 12.5])
set(gca,'Xtick',[0 1 2 3])
set(gca,'Xticklabel',{})

%% Step 4: Calculate dSNP/isolate for each time point, other than the first one
% Test for TP5
RP1=100;RP2=50;
ddSNP=[];
dSNP={};
for ltp = 2:snum;      % Use time points > 1
    figure(1045);hold on;
    set(gca,'LineWidth',1);
    set( gca, 'TickDir', 'out' );
    ylabel('dSNP/isolate');
    SubSetB = find(TPS{ltp}>0);
    SubSetA = find(TP1>0);          % SubsetA is what to compare with
    antnts = ancnti(goodpos);
    dSNP={};        % delta SNP
    mdSNP=[];
    
    % Step 4.2
    % Basis, to make x->1, x>1, x is later time point
    antnts = ancnti(goodpos);
    dSNP={};mdSNP=[];
    for basis = length(SubSetA);
        dSNP{ltp}=[];
        % 1. Pick a subset of TPB; 2. Find how many new SNPs are there when
        % adding a new isolate
        ltp
        basis
        
        BasisIso = SubSetA;
        Remains = SubSetB;

        BasisNT = majorNT(goodpos,BasisIso);
        for rpt2 = 1:length(Remains);
            tmpdSNP = 0;
            rdmRemains = Remains(rpt2);
            newNT = majorNT(goodpos,rdmRemains);
            for pos = 1:length(goodpos);
                if length(unique(BasisNT(pos,:)))==1;
                    if unique(BasisNT(pos,:)) ~= newNT(pos) & newNT(pos) ~=antnts(pos);
                        tmpdSNP = tmpdSNP + 1;
                    end
                end
            end
            dSNP{ltp} = [dSNP{ltp}; tmpdSNP];
        end
    end
    ddSNP=[ddSNP mean(dSNP{ltp})]; 
end

%% Calculate molecular clock
figure(666);hold on;
DIF=ddSNP;
% DIF=ddSNP(2,:);

for ltp = 2:snum;
    h=plot(t(ltp)/365,DIF(ltp-1),'o','Color',c(ltp,:),'MarkerSize',9);set(h,'MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor',c(ltp,:));
end

%%
ty = t/365;
%[S Srange c d e] = regress(DIF' ,[ty(2:length(ty)) [1;1;1;1;1]])

[B,S] = polyfit(ty(2:length(ty)), DIF', 1);

plot([0 2],[0*B(1)+B(2) 2*B(1)+B(2)],'-', 'Color', [.75 .75 .75], 'Linewidth',3);
% plot([0 2],[0*Srange(1)+DIF(1) 2*Srange(1)+DIF(1)],'--', 'Color',(c1+c2+c3)/3, 'Linewidth',2);
% plot([0 2],[0*Srange(2)+DIF(1) 2*Srange(2)+DIF(1)],'--', 'Color',(c1+c2+c3)/3, 'Linewidth',2);

xlim([0,2])
ylim([0 2])
% set(gca,'YTick', [0:.2:.8]);
% set(gca,'XTick', [0.8 , 2.25]);
% set(gca,'XTickLabel',{'',''});
set(gca,'Fontsize',15)
set( gca, 'TickDir', 'out' );
set(gca,'linewidth',1)
end



