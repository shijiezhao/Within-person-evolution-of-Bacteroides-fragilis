function p = generate_positions(SampleDirs, SampleNames,maxFQ, parallel, jobsubmitoptions,TEMPORARYFOLDER,SCRIPTSDIRECTORY)



load for_matlab_scripts
[ChrStarts,GenomeLength,~,~]=genomestats(REF_GENOME_DIRECTORY);


timesvariant=zeros(GenomeLength,1);



%run analysis
cmds={};
for i=1:length(SampleDirs)
    if parallel == 1
        cmds{end+1}=['matlab -r "path(' char(39) SCRIPTSDIRECTORY char(39) ',path); generate_positions_single_sample(' char(39) SampleDirs{i} char(39) ',' char(39) SampleNames{i} char(39) ',' num2str(maxFQ) ');"'];
    else
        generate_positions_single_sample(SampleDirs{i},SampleNames{i},maxFQ);
    end
end

send_jobs_to_cluster(cmds,jobsubmitoptions,parallel,{'.'});



%load files
for i=1:length(SampleDirs)
    %http://www.vsoch.com/2010/11/loading-dynamic-variables-in-a-static-workspace-in-matlab/
    pos=load([TEMPORARYFOLDER '/vcf_' SampleNames{i} '.mat']);
    if numel(pos.Positions)>2 % HC 9/13/2013
        x=chrpos2index(pos.Positions,ChrStarts);
        timesvariant(x)=timesvariant(x)+1;
    end
    delete([TEMPORARYFOLDER '/vcf_' SampleNames{i} '.mat'])
end

p=find(timesvariant>0 & timesvariant <numel(SampleDirs));

fprintf(['Not considering ' num2str(sum(timesvariant==numel(SampleDirs))) ' positions where all samples have a variant compared to the reference...\n'])


