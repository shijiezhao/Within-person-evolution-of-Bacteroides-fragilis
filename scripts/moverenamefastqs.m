function cmds = moverenamefastqs(IsolateTable)

%Modified by TDL on 10/14/2013. Now the first column in samples.csv ('Batch') must
%contain either (1) entire path to the folder containing all fastq files or
% (2) entire path to folder containing subfolder (1 subfolder per sample)

% Modified by TDL on 10/14/2014 to accomdate the same sample being split
% across multiple batches. Requires paired end reads


cmds = {} ;

k=0;
forwardcmd=[];
reversecmd=[];

for i=1:numel(IsolateTable)
    
    s = IsolateTable(i) ;
    if ~exist(s.Sample,'dir')
        mkdir(s.Sample)
    end
    
    
    if (~exist([s.Sample '/' s.Sample '_1.fastq'],'file') | ~exist([s.Sample '/' s.Sample '_2.fastq'],'file')) ...
            & ~(exist([s.Sample '/sickle2050/pairedendtb/'],'dir') | exist([s.Sample '/sickle2050/paired_ecoli_tet_simple/'],'dir') | exist([s.Sample '/sickle2050/pairedendecoli/'],'dir'));
        
        k=k+1;
        sourcefolders = s.Batch;
        
        
        forwardcmd=[forwardcmd ' cat '];
        reversecmd=[reversecmd ' cat '];
        
        
        %add instructions to copy each raw file
        for j=1:numel(sourcefolders)
            
            %process folder name
            sourcefolder=sourcefolders{j};
            if sourcefolder(end)=='/'
                sourcefolder=sourcefolder(1:end-1);
            end
            
            %find fastqs
            if exist([sourcefolder '/' s.ProviderName],'dir')
                forward=dir([sourcefolder '/' s.ProviderName '/*1.fastq']); % HC fixed bug
                reverse=dir([sourcefolder '/' s.ProviderName '/*2.fastq']); % HC fixed bug
                sourcefolder=[sourcefolder '/' s.ProviderName];
            else
                forward=dir([sourcefolder '/' s.ProviderName '*1.fastq']); % TL fixed back
                reverse=dir([sourcefolder '/' s.ProviderName '*2.fastq']); % TL fixed back
            end
            
            %check
            if (numel(forward) ~= 1 | numel(reverse) ~= 1)  & ~exist([s.Sample '/' s.Sample '_' num2str(1) '.fastq'],'file')
                error(['Could not find exactly 1 fastq file in each direction for: ' s.ProviderName  ' in batch:' sourcefolder]);
            end
            
            forwardcmd = [forwardcmd ' ' sourcefolder '/' forward.name];
            reversecmd = [reversecmd ' ' sourcefolder '/' reverse.name];
            
        end
        
        if k==100
            cmds{end+1}=[forwardcmd ' > ' s.Sample '/' s.Sample '_1.fastq'];
            cmds{end+1}=[reversecmd ' > ' s.Sample '/' s.Sample '_2.fastq'];
            k=0;
            forwardcmd=[];
            reversecmd=[];
        else
            forwardcmd=[forwardcmd ' > ' s.Sample '/' s.Sample '_1.fastq ; '];
            reversecmd=[reversecmd ' > ' s.Sample '/' s.Sample '_2.fastq ; '];
        end
        
    end
    
    
end

if k~=0
    cmds{end+1}=forwardcmd;
    cmds{end+1}=reversecmd;
end
