% eliminate homologous regions between MEDs within L01
if strcmp(Donor,'S01'); chrtorm=[36,44,46,89,142,164,348,61];end;
%% 3. Load files for processing
cd(workingdir)
load('candidate_mutation_table')
load('coveragematrix')

%% 4. Remove undesired samples based on name and/or coverage
% 4.1. Coverage calculation: counts(39:candidats:samples)
coverage=squeeze(sum(counts(1:8,:,:)));             % Total counts
coverage_forward=squeeze(sum(counts(1:4,:,:)));     % Forward counts
coverage_reverse=squeeze(sum(counts(5:8,:,:)));     % Reverse counts
% 4.2. Filter samples based on coverage
goodsamples = mean(coverage) > min_average_coverage_to_include_sample;
goodsamples(rmsamples) = 0;                         % Remove bad samples
goodsamples = goodsamples>0;

%% Update sample names, blabla
AllSampNames=SampleNames;
SampleNames=SampleNames(goodsamples);
Allcounts=counts;
counts=counts(:,:,goodsamples);
AllQuals=Quals;
Quals=Quals(:,goodsamples);
Allcoverage=coverage;
coverage=coverage(:,goodsamples);
Nsample=numel(SampleNames);
Allcoverage_forward=coverage_forward;
coverage_forward=coverage_forward(:,goodsamples);
Allcoverage_reverse=coverage_reverse;
coverage_reverse=coverage_reverse(:,goodsamples);
Allin_outgroup=in_outgroup;
in_outgroup=in_outgroup(goodsamples);
num_in_outgroup = sum(in_outgroup);
Quals = -1*Quals; %use -Quals because this way higher numbers are more confident
Allall_coverage_per_bp = all_coverage_per_bp;
all_coverage_per_bp = all_coverage_per_bp(goodsamples,:);

%% 5. Extract the reference genome information
% 5.1. Chromosome information, genome length, Scaffold names
[ChrStarts, GenomeLength, ~, ScafNames]= genomestats(REFGENOMEFOLDER);
fr = fastaread([REFGENOMEFOLDER '/genome.fasta']);
% 5.2. Extract reference genotype at the candidate locations
refnt = extract_outgroup_mutation_positions(REFGENOMEFOLDER, p2chrpos(p,ChrStarts));
% 5.3. Make some basic structures for finding mutations
[majorAF, majorNT, minorNT, minorAF] = div_major_allele_freq(counts);
[~,refnti]=ismember(refnt,NTs);

%% 6. Find positions with fixed mutations
contig_positions=p2chrpos(p, ChrStarts);            % Contig and position of each candidate SNP
% 6.1. Find the nucleotide identity at each position.
Calls=majorNT;
% 6.2. Delete sample*candidates based on 
Calls(Quals < min_qual_for_call)=0;                 % 1). quality score
Calls(majorAF< min_maf_for_call)=0;                 % 2). Major allele frequency
Calls(coverage_forward < min_cov_for_call_per_strand)=0;    % 3). forward coverage
Calls(coverage_reverse < min_cov_for_call_per_strand)=0;    % 4). reversed coverage
Calls(sum(Calls(:,~in_outgroup)<1,2)>=(sum(~in_outgroup)*max_fraction_ambigious_samples),:)=0;  % 5). Candidates with too many ambigous positions
Calls(median(coverage(:,~in_outgroup),2)<min_median_coverage_position,:)=0;                     % 6). Coverage

% 6.2. Determine the major allele in the outgroup

% First round: use closest 
for i = 1:size(SampleNames);
    if strcmp(closest,SampleNames(i));
        ancnti=Calls(:,i);
        fprintf 'The closest strain is', SampleNames(i)
    end
end


%% Read and stitch the mobile element sequences
read_mb_region()


%% 6.3. ??????
% Mutual quality analysis
% mutatedp = find_mutated_positions(Calls(:,~in_outgroup)) ;  
[mutatedp, MutQualIsolates] = ana_mutation_quality(Calls(:,~in_outgroup),Quals(:,~in_outgroup))
% Filter based on Mutual quality ?????????
hasmutation=((Calls~=repmat(refnti,1,Nsample)) & Calls>0 & repmat(mutatedp,1,Nsample)>=FQ_cutoff);
chrs = p2chrpos(p,ChrStarts);
goodpos=find(sum(hasmutation,2)>0 & ~(ismember(chrs(:,1),chrtorm)));

%% Remove SNPs in the mobile regions
trm=[];
for i = 1:size(newstartend,1);
    A = find(p(goodpos)>newstartend(i,1) & p(goodpos)<newstartend(i,2));
    trm=[trm;A];
end
goodpos = setdiff(goodpos,goodpos(trm));

% goodsamples = sum((Calls(goodpos,:)>0),1)> (1-max_ambigious_locations) * sum(numel(goodpos));

%% Second round!

% Rules: (1). Parsimony; (2). Majority from the outgroups; (3). Minimize
% variation

for i = 1:size(SampleNames);
    if strcmp(closest,SampleNames(i));
        ancnti=Calls(:,i);
        fprintf 'The closest strain is', SampleNames(i)
    end
end


if step == 2;
    for i = 1:size(SampleNames);
        if strcmp(closest,SampleNames(i));
            ancnti=Calls(:,i);
            fprintf 'The closest strain is', SampleNames(i)
        end
    end
    for i = 1:length(goodpos);
        if dnapar_ant(i) > 0;                   % Use the parsimonious nt
            ancnti(goodpos(i)) = dnapar_ant(i);
        end
        if ancnti(goodpos(i)) == 0;             % Still, use the consensus of the outgroups
            tmp = Calls(goodpos(i),in_outgroup>0);
            if length(tmp(tmp>0))>0;
                ancnti(goodpos(i))=mode(tmp(tmp>0));
            end
        end
        
    end
end

% Fill in ancestral allele with nt that minimize dmrca if it's still zero
if step==2;
Calls_for_treei=majorNT(goodpos,samplestoplot);
Calls_for_treei = [ancnti(goodpos), Calls_for_treei];
    while sum(Calls_for_treei(:,1)>0)<size(Calls_for_treei,1); 
        % When there is unclear position, use the one that minimize variation
        tmp_call_for_treei = Calls_for_treei(Calls_for_treei(:,1)>0,:);
        A = find(Calls_for_treei(:,1)==0); 
        A = A(1);           % Only work on the first unclear position
        nts = unique(Calls_for_treei(A,:)); nts = nts(nts>0);
        nt1 = nts(1); nt2 = nts(2);
        tmp1 = [Calls_for_treei(A,:);tmp_call_for_treei]; tmp1(1,1)=nt1;
        tmp2 = [Calls_for_treei(A,:);tmp_call_for_treei]; tmp2(1,1)=nt2;
        dMRCA_tmp1 = []; dMRCA_tmp2 = [];
        for i = 2:size(tmp1,2);
            dMRCA_tmp1 = [dMRCA_tmp1; sum(tmp1(:,1)~=tmp1(:,i))];
            dMRCA_tmp2 = [dMRCA_tmp2; sum(tmp2(:,1)~=tmp2(:,i))];
        end
        if var(dMRCA_tmp1)<var(dMRCA_tmp2);
            ancnti(goodpos(A))=nt1;
        else;
            ancnti(goodpos(A))=nt2;
        end
        Calls_for_treei=majorNT(goodpos,samplestoplot);
        Calls_for_treei = [ancnti(goodpos), Calls_for_treei];
    end

hasmutation=((Calls~=repmat(ancnti,1,Nsample)) & Calls>0 & repmat(mutatedp,1,Nsample)>=FQ_cutoff);

end

%% Move the files to a folder, so that the working directory is not crowded
% Move all previous tree files to a storage directory 'PhylogeneticTrees'
thisyear=year(datetime('today'));
eval(['! mv ' num2str(thisyear) '* PhylogeneticTrees'])
eval('! rm outtree'); eval('! rm outfile'); eval('! rm temp.sh');
% Use dnapar to construct a parsimony tree.
plot_phylogeny();
