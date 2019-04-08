%% 0. Purposes of this script
% 0.1. Identify de novo mutations
% 0.2. Find large in-dels (potential HGTs)
% 0.3. Calculate dN/dS ratio for the whole sample
% 0.4. Visualization of the results

%% 1. Initial Parameters
% Much of this will need to vary according to your reference genome,
% coverage, particular samples, and how well the experiment was done.
% 1.1. Parameters for finding de novo mutations between samples
min_average_coverage_to_include_sample = 5;      % Filter out samples have lower than that number of coverage
max_fraction_ambigious_samples = .90;            % If more than x% of the samples have ambiguous NT, discard the candidate location
min_median_coverage_position = 10;               % Remove candidate locations have lower than this coverage
min_qual_for_call = 10;                          % Remove sample*candidate that has lower than this quality
min_maf_for_call = .6;                           % Remove sample*candidate
min_cov_for_call_per_strand = 1;                 % Remove sample*candidate
FQ_cutoff=60;      %min -FQ that samples supporting both alleles must have
% 1.3. Promoter parameter: how far upstream of the nearest gene to annotate
promotersize=300;

%% 2. Enviornment set up -- probably won't need to change
% 2.1. Set up directories
workingdir=char(pwd);   % set up current directory as working dir.
cd('../')                % move up one directory
masterdir=char(pwd);    % set up the parent dir. as masterdir
REFGENOMEFOLDER=[masterdir '/references/BfragCR'];  
SCRIPTSDIRECTORY = [masterdir '/scripts'];
path(SCRIPTSDIRECTORY,path);    % Add the srcipt folder to the search path
TreefileName = 'BfragCR';       % Tree name
SavefileName = 'TEST.mat';
% 2.2. Initilize ATCG
NTs='ATCG';

%% 3. Load files for processing
cd(workingdir)
cd('public_data_cross_lineage')
load('candidate_mutation_table')
load([REFGENOMEFOLDER '/cds_sorted.mat'])


%% 4. Remove undesired samples based on name and/or coverage
% 4.1. Coverage calculation: counts(39:candidats:samples)
%%%%%%%%%%%%%::: 39?
coverage=squeeze(sum(counts(1:8,:,:)));             % Total counts
coverage_forward=squeeze(sum(counts(1:4,:,:)));     % Forward counts
coverage_reverse=squeeze(sum(counts(5:8,:,:)));     % Reverse counts
% 4.2. Filter samples based on coverage
goodsamples = mean(coverage) > min_average_coverage_to_include_sample;
SampleNames=SampleNames(goodsamples);
counts=counts(:,:,goodsamples);
Quals=Quals(:,goodsamples);
coverage=coverage(:,goodsamples);
Nsample=numel(SampleNames);
in_outgroup=in_outgroup(goodsamples);
coverage_forward=coverage_forward(:,goodsamples);
coverage_reverse=coverage_reverse(:,goodsamples);
Quals = -1*Quals; %use -Quals because this way higher numbers are more confident

%% 5. Extract the reference genome information
% 5.1. Chromosome information, genome length, Scaffold names
[ChrStarts, GenomeLength, ~, ScafNames]= genomestats(REFGENOMEFOLDER);
% 5.2. Extract reference genotype at the candidate locations
refnt = extract_outgroup_mutation_positions(REFGENOMEFOLDER, p2chrpos(p,ChrStarts));
% 5.3. Make some basic structures for finding mutations
[majorAF, majorNT, minorNT, minorAF] = div_major_allele_freq(counts);
[~,refnti]=ismember(refnt,NTs);


%% 6. Find positions with fixed mutations
% 6.1. Find the nucleotide identity at each position.
Calls=majorNT;
% 6.2. Delete sample*candidates based on 
% 1). quality score
% 2). Major allele frequency
% 3,4). forward/reversed coverage
% 5). Candidates with too many ambigous positions
% 6). Coverage
Calls(Quals < min_qual_for_call)=0;
Calls(majorAF< min_maf_for_call)=0;
Calls(coverage_forward < min_cov_for_call_per_strand)=0;
Calls(coverage_reverse < min_cov_for_call_per_strand)=0;
Calls(sum(Calls(:,~in_outgroup)<1,2)>=(sum(~in_outgroup)*max_fraction_ambigious_samples),:)=0;
Calls(median(coverage(:,~in_outgroup),2)<min_median_coverage_position,:)=0;
% 6.2. Determine the major allele in the outgroup
ancnti=mode(Calls(:,in_outgroup>0),2);
ancnti(ancnti==0)=refnti(ancnti==0); % if empty, use the reference
ancnti=refnti;
% 6.3. ??????
% Mutual quality analysis
[MutQual, MutQualIsolates] = ana_mutation_quality(Calls(:,~in_outgroup),Quals(:,~in_outgroup)) ;  
% Filter based on Mutual quality ?????????
fixedmutation=((Calls~=repmat(ancnti,1,Nsample)) & Calls>0 & repmat(MutQual,1,Nsample)>=FQ_cutoff);
hasmutation= fixedmutation; %| diversemutation;
tmpgoodpos=find(sum(hasmutation,2)>0);
%% Only consider variable nucleotides that are from the selected genes
subGL=[];
GL={'BF0864','BF0893','BF1802','BF1803','BF2942','BF3581','BF4056','BF0188','BF1708','BF2848','BF0991','BF1174','BF2755','BF3560'}
% GL={'BF4056'}
for i = 1:length(CDS{1});
    for j = 1:length(GL);
        if strcmp(CDS{1}(i).oldlocustag,GL{j})==1;
            subGL=[subGL;[CDS{1}(i).loc1,CDS{1}(i).loc2]];
        end
    end
end
% Get only information of SNPs from these genes
goodpos=[];
isgoodpos=[];
for i = 1:length(tmpgoodpos);
    if sum((p(tmpgoodpos(i))>subGL(:,1)).*(p(tmpgoodpos(i))<subGL(:,2))) > 0;
        goodpos=[goodpos;tmpgoodpos(i)];
    end
end

%% 7. Create annotation file
positions_to_show_in_table=goodpos;
annotations = annotate_mutations_gb(p2chrpos(p(positions_to_show_in_table),ChrStarts),REFGENOMEFOLDER) ;
annotation_full= append_annotations(annotations, ancnti(positions_to_show_in_table), Calls(positions_to_show_in_table,~in_outgroup), counts(:,positions_to_show_in_table,~in_outgroup), hasmutation(positions_to_show_in_table,~in_outgroup), promotersize) ; %% adds information about particular mutations observed, based on

%% Compute dNdS and CI for each gene
TEMP=[];
for i = 1:length(GL);   % For each gene on the list
    tmp = []; N=0;S=0;  % Record 
    for j = 1:length(annotation_full);
        if strcmp(annotation_full(j).oldlocustag,GL{i})==1;
            tmp = [tmp;j];
            if annotation_full(j).type=='N'; N=N+1;end;
            if annotation_full(j).type=='S'; S=S+1;end;
        end
    end
    subANNOT = annotation_full(tmp);    % the annotation information of SNPs from this gene
    eNS = expectedNS(subANNOT);         % expected dNdS
    [XX,YY,ZZ]=binomialCIdNdS(N,S,eNS);
    TEMP=[TEMP;[i,ZZ,XX,YY]];
end
%% Plot the figure
figure(1);hold on;xlim([0 14]);
Grey1=[64 64 64]/256;G2=[0.8 0.8 0.8]; c1 = [193 220 216]/256;c3 = [212 160 156]/256;c2 = [50 92 127]/256;c4 = [234 201 157]/256;c5 = [139 112 165]/256;c6 = [64 64 64]/256;c7 = (c1+c2+c3)/3;
% Plot bars
for i=1:6;
    h=bar(i,log(TEMP(i,2)),0.5);set(h,'FaceColor',G2);set(h,'LineWidth',0.05);set(h,'EdgeAlpha',0);h(1).BaseValue = log(1);
end
% do not plot the mutations from BF4056, because no inter-lineage mutations
% were found
for i=8:14;
    h=bar(i-1,log(TEMP(i,2)),0.5);set(h,'FaceColor',G2);set(h,'LineWidth',0.05);set(h,'EdgeAlpha',0);h(1).BaseValue = log(1);
end
% Plot baseline for inter-lineage SNPs
plot([0 17],[log(0.1629) log(0.1629)],'--','Color',c2,'LineWidth',1);plot([0 17],[log(1) log(1)],'-','Color',c2,'LineWidth',1)
% Plot the CI
for i=1:6;
    plot([i i], [max(-4,log(TEMP(i,3))),max(-4,log(TEMP(i,4)))],'-','Color',Grey1,'LineWidth',1);
end
for i=8:14;
    plot([i-1 i-1], [max(-4,log(TEMP(i,3))),max(-4,log(TEMP(i,4)))],'-','Color',Grey1,'LineWidth',1);
end

% Format
set(gca,'Fontsize',15);set( gca, 'TickDir', 'out' );;set(gca,'linewidth',1)
set(gca,'Ytick',[log(0.06) log(0.25) log(1) log(4)]);set(gca,'Yticklabel',{'0.06' '0.25' '1' '4'})
ylim([log(0.025) log(4)]);set(gca,'Xtick',1:length(GL));set(gca,'Xticklabel',{})

