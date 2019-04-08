function indelp=score_all_indels(SampleDirs,SampleNames,RefGenome, jobsubmitoptions)

%Tami Lieberman 2013


eval('! rm merged.dindel_output.variants.txt');


global RUN_ON_CLUSTER SCRIPTSPATH CLUSTERDIR;
if RUN_ON_CLUSTER == 1
    mainfolder=CLUSTERDIR;
else
    mainfolder='/Volumes/sysbio/KISHONY LAB/illumina_pipeline';
end

scafdict = makescafdict(RefGenome);
[ChrStarts,GenomeLength,ChromosomeIndicator]=genomestats(RefGenome, 1);
NSample=numel(SampleDirs);

%Combine variant lists
!touch merged.dindel_output.variants.txt
for i=1:NSample
    
    %First filter
    eval(['!python ' SCRIPTSPATH '/dindel-1.01-python/selectCandidates.py'...
        ' -i ' SampleDirs{i} '/sample.dindel_output.variants.txt -o ' SampleDirs{i} ...
        '/sample.mincount.variants.txt --minCount 5']);
    
    %Then combine
    eval(['! cat ' SampleDirs{i} '/sample.mincount.variants.txt >> merged.dindel_output.variants.txt']);
    
end


fprintf(1,'Sorting variants...\n');
%Sort variant list
%Done in this way to allow variables in the executable line
eval(['!/opt/dindel-1.01/dindel --analysis realignCandidates --varFile merged.dindel_output.variants.txt '...
    '--outputFile sorted.mincount --ref ' mainfolder '/Reference_Genomes/' RefGenome '/genome.fasta']);
%Sort list into windows
%Done in this way to allow variables in the executable line

if exist('indelwindows', 'dir')
    eval(['!rm -r indelwindows']);
end
mkdir('indelwindows');
eval(['!python ' SCRIPTSPATH '/dindel-1.01-python/makeWindows.py --inputVarFile sorted.mincount.variants.txt --windowFilePrefix indelwindows/window --numWindowsPerFile 1000']);

%Count number of windows
numwindows=numel(dir('indelwindows/window*'));

fprintf(1,'Realigning in variant windows...\n');

%Run each window on each strain
cmds={}; j=0; cmd=[];
for i=1:NSample
    
    eval(['!touch ' SampleNames{i} '/sample.dindel_stage2_outputfiles.txt']);
    
    for w=1:numwindows
        
        cmd=[cmd '/opt/dindel-1.01/dindel --analysis indels --doPooled  --bamFile '...
            SampleDirs{i} '/aligned.sorted.bam --ref  ' CLUSTERDIR '/Reference_Genomes/'...
            RefGenome '/genome.fasta --varFile indelwindows/window.' num2str(w) '.txt'...
            ' --libFile ' SampleDirs{i} '/sample.dindel_output.libraries.txt'...
            ' --outputFile ' SampleNames{i} '/sample.dindel_stage2_output_windows_pooled.' num2str(w) ' --quiet;'];
        
        j=j+1;
        if j> 5
            cmds{end+1}=cmd;
            j=0;
            cmd=[];
        end
        
        eval(['!echo ' SampleNames{i} '/sample.dindel_stage2_output_windows_pooled.' ...
            num2str(w) '.glf.txt >> ' SampleNames{i} '/sample.dindel_stage2_outputfiles.txt']);
    end
    
    
end

cmds{end+1}=cmd; %% add last one

run_parallel_unix_commands_fast(cmds,'short -W 2:00',1,{pwd});


fprintf(1,'Processing indel results...\n');

%merge for each sample
cmds={};
for i=1:NSample
    cmds{end+1}=['python ' SCRIPTSPATH '/dindel-1.01-python/mergeOutputPooled.py '...
        '--inputFiles ' SampleNames{i} '/sample.dindel_stage2_outputfiles.txt '...
        '--outputFile ' SampleNames{i} '/pooledvariantCalls.VCF --ref ' mainfolder '/Reference_Genomes/'...
        RefGenome '/genome.fasta --numSamples 1 --numBamFiles 1'];
end

run_parallel_unix_commands_fast(cmds,'mini -W 10',1,{pwd});


fprintf(1,'Combining indel results...\n');
fprintf(1,'If this part runs slowly, need to parallelize the end of score_all_indels...\n');

%combine to make list of ps for returning
potentialp=zeros(100000,1);

c=1;
for i=1:NSample
    
    fprintf(1,[SampleNames{i} '...']);
    vcf = read_vcf_file_dindel([SampleNames{i} '/pooledvariantCalls.VCF']) ;
    if isstruct(vcf)
        % vcf= vcf([vcf.nf]>1 & [vcf.nr]>1);
        indel_lengths=max([cellfun(@numel,{vcf.ref}); cellfun(@numel,{vcf.alt})]);
        
        for j=1:numel(vcf)
            chr=scafdict(vcf(j).scaf);
            pos=vcf(j).pos;
            start=double(ChrStarts(chr)+pos);
            potentialp(c:c+indel_lengths(j))=start:start+indel_lengths(j);
            c=c+indel_lengths(j);
        end
    end
end


indelp=unique(potentialp);





