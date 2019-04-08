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
cd('All_lineages')
load('candidate_mutation_table.mat')


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
goodpos=find(sum(hasmutation,2)>0);



%% Prepare files for tree building
samplestoplot=[1:numel(SampleNames)];
quality_positions = (MutQual>FQ_cutoff & sum(Calls==0,2)<(.2*numel(SampleNames)) & median(coverage,2)>5);
quality_positions = goodpos;
Calls_for_treei=Calls(quality_positions,samplestoplot);
calls_for_tree=zeros(size(Calls_for_treei));
calls_for_tree(Calls_for_treei>0)=NTs(Calls_for_treei(Calls_for_treei>0));
calls_for_tree(Calls_for_treei==0)='N';
% ADD REFERENCE AT THESE POSITIONS
outgroup_nts = extract_outgroup_mutation_positions(REFGENOMEFOLDER, p2chrpos(p(quality_positions),ChrStarts));
TreeSampleNames= {SampleNames{samplestoplot}};
calls_for_tree = [outgroup_nts, calls_for_tree];
TreeSampleNames={'Reference' SampleNames{samplestoplot}};



%% Draw NJ tree
Seqs = {}; sel=[]; % Temporarily store positios that are selected
Names={};
tmpi = Calls_for_treei;
tmp = calls_for_tree;

% Delete bad samples
sel=[];
for i = 1:size(tmpi,2); %% Delete samples with too many 0's
    badalel = sum(tmpi(:,i)==0);
    if badalel < 10000;  % delete samples with more than 10000 positions that have zero's
        sel = [sel i];
    end
end

k=[1 sel+1];  % translate the selected samples into calls_for_tree matched
Names={TreeSampleNames{k}}; % select sample names
tmpi = tmpi(:,sel); % delete samples with more than 10000 0's
tmp = tmp(:,k); % delete from the calls_for_tree file


% Delete positions with too many 0's
sel=[]; 
for i = 1:size(tmpi,1);
    if sum(tmpi(i,:)>0) > size(tmpi,2)-5;  %
        sel=[sel;i];
    end
end
% sel=find(prod(Calls_for_treei,2)>0);
tmp = tmp(sel,:);
tmpi = tmpi(sel,:);


% Delete bad samples
sel=[];
for i = 1:size(tmpi,2); %% Delete samples with too many 0's
    badalel = sum(tmpi(:,i)==0);
    if badalel < 2000;  % delete samples with more than 10000 positions that have zero's
        sel = [sel i];
    end
end

k=[1 sel+1];  % translate the selected samples into calls_for_tree matched
Names={Names{k}}; % select sample names
tmpi = tmpi(:,sel); % delete samples with more than 10000 0's
tmp = tmp(:,k); % delete from the calls_for_tree file


%% N-J computations
% for sample_number = 1:size(tmp,2);
%     Seqs{sample_number} = char(tmp(:,sample_number));
% end
D=[];
tmpi_ci = tmp(:,1);
tmpi_ci = 1*(tmpi_ci==65)+2*(tmpi_ci==84)+3*(tmpi_ci==67)+4*(tmpi_ci==71);
tmpii = [tmpi_ci tmpi];

for i = 1:size(tmpii,2)-1;
    i
    for j = i+1:size(tmpii,2);
        Seqs={NTs(ones(size(tmpii,1),1)) NTs(ones(size(tmpii,1),1))};
        diff = find(tmpii(:,i)>0 & tmpii(:,j)>0);
        Seqs{1}(diff) = NTs(tmpii(diff,i));
        Seqs{2}(diff) = NTs(tmpii(diff,j));
        D=[D seqpdist(Seqs)];

    end
end

% N-J
PhyloTree = seqneighjoin(D,'equivar',Names)
% Save the tree
phytreewrite('Figures/Figure_1a.tree', PhyloTree,'BranchNames',0)

%% Plot supplementary figure 1b

DIF=[];

NCTC9343_Intra_distance={[],[],[],[],[],[],[],[],[],[],[],[]}
for i = 2:size(tmpi,2)-1;
    for j = i+1:size(tmpi,2);
        diff = sum(tmpi(:,i)>0 & tmpi(:,i)~=tmpi(:,j) & tmpi(:,j)>0);
        if diff<1000;
           ID=strsplit(Names{i},'_');
           ID = ID{1};
           subjectid = str2num(ID(2:3));
           NCTC9343_Intra_distance{subjectid}=[NCTC9343_Intra_distance{subjectid} diff];
        end
        %         Bugs=[Bugs;diff];
%         if strcmp(SampleNames{i}(1:3),'S02')==0;
            DIF=[DIF;diff];
%             SampleNames{i}
%             SampleNames{j}
%         end
    end
end

%% Plot the histogram: sup_figure 1b
c1=[175 65 27]/255;c2=[139 175 198]/255;
figure(1); hold on;
inter_dif = DIF(DIF>500);   
intra_dif = DIF(DIF<500);% This is empirical, but true. Pairs with distance <500 are from a same lineage
% Plot the histogram
pts=0000:1000:25000;
[bincounts] = histc(inter_dif,pts);
h=bar(pts,bincounts,'histc');
h.FaceColor = c2;

[bincounts] = histc(intra_dif,pts);
h=bar(pts,bincounts,'histc');
h.FaceColor = c1;

int_57 = ones(29,1)*18500;
[bincounts] = histc(int_57,pts);
h=bar(pts,bincounts,'histc')
h.FaceColor = c1;
set( gca, 'TickDir', 'out' );
box off

xlim([-500 20000])
