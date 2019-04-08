function [p,coveragethresholds] = find_diverse_positions_no_control(params, SampleDirs, SampleNames, parallel, jobsubmitoptions, TEMPORARYFOLDER,SCRIPTSDIRECTORY)


%% Tami Lieberman 2012-2016

%p is a structure that contains the list of genomic positions meetings
%params

%coveragethresholds -- contains cdf cutoffs .01:.01:1.0.
%only counts positions where at least 1 read aligned
%purpose of this data structure is to allow downstream processes to
%remove positions with excess coverage in a way that accounts for
%variation in coverage between samples


%%



load for_matlab_scripts
[~,GenomeLength,~,~]=genomestats(REF_GENOME_DIRECTORY);


if ~exist('find_diverse_positions_log')
    mkdir('find_diverse_positions_log')
end

% create a dummy sample because find_diverse_positions_single_sample is written to include an isogenic control 
MAF_control=[0];
save('for_finding_diverse_positions', 'params', 'GenomeLength', 'MAF_control', '-v7.3')

% run each non-control sample
cmds={};
for i=1:size(SampleNames)
    if parallel
        cmds{end+1}=['matlab -r "path(' char(39) SCRIPTSDIRECTORY char(39) ',path); find_diverse_positions_single_sample(' char(39) SampleDirs{i} char(39) ',' char(39) SampleNames{i} char(39) ',' char(39)  TEMPORARYFOLDER char(39) ');"'];
    else
        find_diverse_positions_single_sample(SampleDirs{i},SampleNames{i}, TEMPORARYFOLDER);
        
    end
end

send_jobs_to_cluster(cmds,jobsubmitoptions,parallel,{'.'});



%load files
p=zeros(GenomeLength,1);
coveragethresholds=zeros(100,numel(SampleNames));

for i=1:size(SampleNames)
    %http://www.vsoch.com/2010/11/loading-dynamic-variables-in-a-static-workspace-in-matlab/
    diverse=load([TEMPORARYFOLDER '/diverse_' SampleNames{i} '.mat']);
    p=p+diverse.p_sample;
    coveragethresholds(:,i)=diverse.coveragethresholds_sample;
    delete([TEMPORARYFOLDER '/diverse_' SampleNames{i} '.mat'])
end

%return all positions potentially diverse in at least one sample
p=find(p>0);



