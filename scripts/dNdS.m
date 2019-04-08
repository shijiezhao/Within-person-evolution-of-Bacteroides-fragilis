
%% 2. Enviornment set up -- probably wont need to change
% 2.1. Set up directories
workingdir=char(pwd);   % set up current directory as working dir.
cd('../../..')                % move up one directory
masterdir=char(pwd);    % set up the parent dir. as masterdir
REFGENOMEFOLDER=[masterdir '/Reference_Genomes/BfragCR'];  
SCRIPTSDIRECTORY = [masterdir '/scripts'];
path(SCRIPTSDIRECTORY,path);
% 2.2. Initilize ATCG
NTs='ATCG';

cd(workingdir)
load('allgenes.across.donors.mat')
%% 9. Find genes mutated more than once
locustags_p={annotation_full(:).locustag};
tags_p=locustags_p(cellfun(@numel,locustags_p)>0);
for i=1:numel(locustags_p)
    if isempty(locustags_p{i})
        locustags_p{i}='intergenic';
    end
end
[number_times_gene_mutated_p,mutated_genes_p]=find_count_duplicates(tags_p);
genes_mutated_multiple_times=mutated_genes_p(number_times_gene_mutated_p>1); %Should consider putting some gene length threshold in here;

genes_mutated_trice=mutated_genes_p(number_times_gene_mutated_p>2);

in_gene_mutated_multiple_times=ismember(locustags_p,genes_mutated_multiple_times);
in_gene_mutated_trice=ismember(locustags_p,genes_mutated_trice);

gene_mutated_in_Person={['BF9343_0949'],['BF9343_1115'],['BF9343_1627'],['BF9343_1720'],['BF9343_0826'],['BF9343_0855'],['BF9343_1721'],['BF9343_3483'],['BF9343_3949'],['BF9343_2669'],['BF9343_2761'],['BF9343_3014'],['BF9343_3467']}
in_gene_mutated_in_person=ismember(locustags_p,gene_mutated_in_Person);

%% 10. dNdS analysis across all mutation
dNdSfilename=[REFGENOMEFOLDER '/dNdStools.mat'];
if exist(dNdSfilename,'file')
    load(dNdSfilename)
else
    div_mutation_type_probability_matrix(REFGENOMEFOLDER, promotersize);
    save(dNdSfilename, 'genomewide_possibilities', 'cds_possibilities')
end
%categories={find(~cellfun(@isempty,{annotation_full.locustag})) find(~cellfun(@isempty,{annotation_full.locustag}) & ~in_gene_mutated_multiple_times) find(in_gene_mutated_multiple_times) find(in_gene_mutated_trice) find(in_gene_mutated_in_person)} ;
categories={find(~cellfun(@isempty,{annotation_full_all.locustag})) find(~cellfun(@isempty,{annotation_full.locustag}) & ~in_gene_mutated_multiple_times) find(in_gene_mutated_multiple_times & ~in_gene_mutated_trice) find(in_gene_mutated_trice) find(in_gene_mutated_in_person)} ;
categories={find(~cellfun(@isempty,{annotation_full_all.locustag})) find(~cellfun(@isempty,{annotation_full.locustag}) & ~in_gene_mutated_trice & ~in_gene_mutated_in_person) find(in_gene_mutated_trice) find(in_gene_mutated_in_person)} ;









category_names={'All genes across donors','Genes mutated once','Genes mutated triple times','Genes mutated in single donors'};

mutationmatrix=zeros(4,4,numel(categories)); %last dimension is context
typecounts=zeros(numel(categories),4); %NSPI %first dimension is context
%%

G_scores_1=[];
for context=1:numel(categories)
    muts_in_category=categories{context};
    
    if context>1;
    for j=1:numel(muts_in_category);
        
        i=muts_in_category(j);
        anc=find(NTs==annotation_full(i).anc);
        new=annotation_full(i).nts;
        
        if context == 2 & annotation_full(i).locustag(1)=='B';
            G_scores_1 = [G_scores_1 1000/sum(ExpectedNumberOfMutationPerGene{str2num(annotation_full(i).locustag(8:11))})];
        end
        if find(new==NTs(anc)) & numel(new)==2
            new(new==NTs(anc))=[];
            new(new=='N')=[];
            new=find('ATCG'==new);
            
            mutationmatrix(anc,new,context)=mutationmatrix(anc,new,context)+1;
            
            if annotation_full(i).type=='N'
                typecounts(context,1)=typecounts(context,1)+1;
            elseif annotation_full(i).type=='S'
                typecounts(context,2)=typecounts(context,2)+1;
            elseif annotation_full(i).type=='P'
                typecounts(context,3)=typecounts(context,3)+1;
            elseif annotation_full(i).type=='I'
                typecounts(context,4)=typecounts(context,4)+1;
            end
        end
    end
    end
    if context==1;
    for j=1:numel(muts_in_category);
        
        i=muts_in_category(j);
        anc=find(NTs==annotation_full_all(i).anc);
        new=annotation_full_all(i).nts;
        
        if find(new==NTs(anc)) & numel(new)==2
            new(new==NTs(anc))=[];
            new(new=='N')=[];
            new=find('ATCG'==new);
            
            mutationmatrix(anc,new,context)=mutationmatrix(anc,new,context)+1;
            
            if annotation_full_all(i).type=='N'
                typecounts(context,1)=typecounts(context,1)+1;
            elseif annotation_full_all(i).type=='S'
                typecounts(context,2)=typecounts(context,2)+1;
            elseif annotation_full_all(i).type=='P'
                typecounts(context,3)=typecounts(context,3)+1;
            elseif annotation_full_all(i).type=='I'
                typecounts(context,4)=typecounts(context,4)+1;
            end
        end
    end
    end
end


mutO = div_matrix2_6types(mutationmatrix);

percentN_matrix=(genomewide_possibilities(:,:,1)./(genomewide_possibilities(:,:,1)+genomewide_possibilities(:,:,2)));
percentN_types=div_matrix2_6types(percentN_matrix)./2;

No=typecounts(:,1);
So=typecounts(:,2);

expectedN=mutO'*percentN_types;
expectedS=sum(mutO)'-expectedN;

NSe=expectedN./expectedS;
NSo=typecounts(:,1)./typecounts(:,2);

%%
%mycolors=([[117	133	148]/256 ; rgb('black');[83	119	122]/256;[192	41	66]/256;[84	36	55]/256;[213/256 94/256 0]; rgb('darkblue'); rgb('red'); rgb('blue'); rgb('black'); rgb('gray')]);

% mycolors=([rgb('grey') ; rgb('black');[51	102	153]/256;[188   56  50]/256;[84	36	55]/256;[213/256 94/256 0]; rgb('darkblue'); rgb('red'); rgb('blue'); rgb('black'); rgb('gray')]);
Red1=[244 113 126]/256; Red2=[212 160 156]/256;
Blue1=[224 236 234]/256; Blue2=[193 219 216]/256; Blue3=[100 142 180]/256;
Grey1=[64 64 64]/256; Grey2=[200 200 200]/256;

c1 = Grey1; c2=Grey2;c3=Blue2;c4=Red2;
mycolors=([c1;c2;c3;c4]);

%Indidiviudal
[upperCI,lowerCI]=binomialCIdNdS_new(No,So,NSe)
% NSo=[nso;NSo];NSe=[nse;NSe];upperCI=[uCI;upperCI];lowerCI=[lCI;lowerCI];
figure(14); clf; hold on;
for i=1:numel(NSo)
    h=bar(i,log(NSo(i)/NSe(i)),0.5);
    set(h,'FaceColor',mycolors(i,:));
%     set(h,'EdgeColor','black');
    set(h,'LineWidth',0.05);
    h(1).BaseValue = log(1);
    set(h,'EdgeAlpha',0);

end
%%
plot([1:numel(NSo);1:numel(NSo)],[log(lowerCI');log(upperCI')],'Color','black','LineWidth',1);
set(gca,'Ytick',[log(0.125) log(0.25) log(.5) log(1) log(2) log(4) log(8) log(16) log(32) log(64) log(128)])
set(gca,'Yticklabel',{'0.125' '0.25', '.5' '1' '2' '4' '8' '16' '32' '64' '128'})
% set(gca,'YScale','log')
% set(gca,'Xtick',1:numel(NSo))
set(gca,'XTick',[])
% set(gca,'FontSize',15)
%set(gca,'Xticklabel', category_names)

binocdf(typecounts(:,1)-1,typecounts(:,1)+typecounts(:,2),expectedN./(expectedN+expectedS),'upper')
ylim([log(0.125) log(128)])
set(gca,'TickDir','out')

saveas(gcf,'dNdS.png')
% pause
set(gca,'XColor','w')
xlim([0.25 4.75])

%%
categories={find(~cellfun(@isempty,{annotation_full_all.locustag})) find(~cellfun(@isempty,{annotation_full.locustag}) & ~in_gene_mutated_multiple_times) find(in_gene_mutated_multiple_times & ~in_gene_mutated_trice) find(in_gene_mutated_trice) find(in_gene_mutated_in_person)} ;
A3=annotation_full(categories{3});
A4=annotation_full(categories{4});
A5=annotation_full(categories{5});

%%
figure(44444)
subplot(2,1,1);hold on;
for i = 0:16;
    freq = sum(expectedNumberThriceMutatedGenes==i)/1000;
    h=bar(i,freq,0.5);
    set(h,'FaceColor',mycolors(2,:));
    set(h,'LineWidth',0.05);
    set(h,'EdgeAlpha',0);

    xlim([-1.25 14])
    ylim([0 1])
    set(gca,'Ytick',[0 0.2 0.4 0.6 0.8 1])
    set(gca,'Xtick',[-1 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
    set(gca,'Xticklabel',{'' '0' '1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' ''})
    set(gca,'TickDir','out')
%     set(gca,'FontSize',15)
end
h=bar(9,1,0.3);
set(h,'FaceColor',mycolors(3,:));
set(h,'EdgeAlpha',0);

for j = -1.23:0.001:-1.01;
    h=bar(j,1,0.001);
    h(1).FaceColor=[1 1 1];
    h(1).EdgeColor=[1 1 1];
end
% set(gca,'Xticklabel',{''  ''})

subplot(2,1,2);hold on;
for i = 0:16;
    freq = sum(expectedNumberGenesMutatedMultipleTimesWithinASubject==i)/1000;
    h=bar(i,freq,0.5);
    set(h,'FaceColor',mycolors(2,:));
    set(h,'LineWidth',0.05);
    set(h,'EdgeAlpha',0);

    xlim([-1.25 14])
    ylim([0 1])
    set(gca,'Ytick',[0 0.2 0.4 0.6 0.8 1])
    set(gca,'Xtick',[-1 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
    set(gca,'Xticklabel',{'' '0' '1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' ''})
    set(gca,'TickDir','out')
end
h=bar(13,1,0.3);
set(h,'FaceColor',mycolors(4,:));
set(h,'EdgeAlpha',0);

for j = -1.23:0.001:-1.01;
    h=bar(j,1,0.001);
    h(1).FaceColor=[1 1 1];
    h(1).EdgeColor=[1 1 1];
end
% set(gca,'Xticklabel',{''  ''})


%%
figure(1234);hold on;
for i = 0:30;
    freq = sum(expectedNumberMultipleMutatedGenes==i)/1000;
    h=bar(i,freq,0.5);
    set(h,'FaceColor',mycolors(2,:));
    set(h,'LineWidth',0.05);
    set(h,'EdgeAlpha',0);

    xlim([-0.25 30.25])
    ylim([0 1])
    set(gca,'Ytick',[0 0.2 0.4 0.6 0.8 1])
    set(gca,'Xtick',[0 5 10 15 20 25 30])
    set(gca,'TickDir','out')
    set(gca,'FontSize',15)
end
h=bar(24,1,0.3);
set(h,'FaceColor',mycolors(1,:));
set(h,'EdgeAlpha',0);


% h.FaceColor=[0.5 0.5 0.5];
% set(h,'Ytick',[0 1 2 3 4 5 6 7 8 9 10 ]);


%% Enrichment for membrane protein mutation
% figure(223);hold on;
[e1,p1]=binofit(74,176);
[e2,p2]=binofit(23,28);

[l1,u1]=binomialCIdNdS_new(74,102,0.789);
[l2,u2]=binomialCIdNdS_new(23,5,0.789);

% h=bar(1,e1,0.3);
% set(h,'FaceColor',[0.5 0.5 0.5]);
% set(h,'FaceAlpha',0.3);
% 
% plot([1,1],[p1(1),p1(2)],'Color','black','LineWidth',1);
% 
% h=bar(2,e2,0.3);
% set(h,'FaceColor',[0.5 0.5 0.5]);
% plot([2,2],[p2(1),p2(2)],'Color','black','LineWidth',1);
% set(gca,'FontSize',15)
% set(h,'FaceAlpha',0.3);
% xlim([0.5 2.5])
% set(gca,'Xtick',[])



