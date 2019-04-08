% Go to the lineage folder
cd(['L' Donor(2:3)]);
% 1.2. Samples to remove
chrtorm=[];rmsamples = [];
%% 2. Set the environment parameters
% 2.1. Set up directories
workingdir=char(pwd);              % set up current directory as working dir.
cd('../..')                        % move up one directory
masterdir=char(pwd);               % set up the parent dir. as masterdir
REFGENOMEFOLDER=[masterdir '/references/' RefGen];  
SCRIPTSDIRECTORY = [masterdir '/scripts'];
path(SCRIPTSDIRECTORY,path);       % Add the srcipt folder to the search path
NTs='ATCG';
% 1.0. Get the general parameters;
SNP_identification_parameters()

%% 3. Define files to save
TreefileName = RefGen;             % Tree name
SavefileName = [Donor '.mat'];     % Save the results

%% 4. Identify de novo SNPs
cd(workingdir)
% 4.1. In the first step, identify variable positions using the paramters,
step = 1;
initial_analysis();
% 4.2. Use a python script to find the confident and ambiguous ancestor
% alleles, from the tree files built with dnapar.
eval(['! python ../../scripts/Find_ancestor_alleles.py']);
load('antnt.mat')
% 4.3. The ambiguous ancestors are determined using method described in the
% paper.
step = 2;
initial_analysis();

%% 5. Make a clickable table to view the identified SNPs. 
% ALso, make annotation_full that contains the infomation of all SNPs
positions_to_show_in_table=goodpos;
annotations = annotate_mutations_gb(p2chrpos(p(positions_to_show_in_table),ChrStarts),REFGENOMEFOLDER) ;
annotation_full= append_annotations(annotations, ancnti(positions_to_show_in_table), Calls(positions_to_show_in_table,~in_outgroup), counts(:,positions_to_show_in_table,~in_outgroup), hasmutation(positions_to_show_in_table,~in_outgroup), promotersize) ; %% adds information about particular mutations observed, based on
order=1:numel(SampleNames);
QualSort=0; %set 1 to show mutations with lowest FQ scores up top, 0 to show in order on the genome
clickable_snp_table(annotation_full, Calls(positions_to_show_in_table,order), counts(:,positions_to_show_in_table,order), SampleNames(order), ScafNames, mutatedp(positions_to_show_in_table), QualSort, all_coverage_per_bp, p(positions_to_show_in_table));
save([Donor 'prokka.mat'],'annotation_full')

%% 6. Get the median coverage of isolates from each sample, for supplementary tables
get_median_coverage_per_time_point();

%% 7: Get the intermediary files for summarizing dMRCA, SNP, d(SNP), tip-2-root overtime
summarize_dmrca_snp_number_dsnp();

%% 8. Save the intermidiary files
nsnp=length(annotation_full);
save('SNP_dMRCA_related_intermediary.mat','dmrca_overtime','t2r_overtime','dMRCA_collector','SNP_collector','dSNP','t','nsnp','Tip2Root');

%dmrca_overtime,stddmrca_overtime: dmarca over time, with mean and std
%dSNP: identified dSNP over time, for individual dots
%t2r_overtime: tip-2-root for individual dots
%t: major time points
%dMRCA_collector,SNP_collector: for collector curves


%% Tree counting plots
tree_counting_sector();


% %% Plot the coverage files of mobile elements & update goodpos
% AnalyzeMobileElements();
% %% Create heatmap for mobile elements presence/absence
% Heatmap();