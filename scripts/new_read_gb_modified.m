function new_read_gb_modified(RefGenomeDirectory)


%% Written by Idan Yelin, Roy Kishony, Tami Lieberman, and Hattie Chung

% May need revisions for your reference genome



%% Intialize


if nargin < 2
    use_old_locus_tag = 0;
end

% Load fasta
fastaFilename = [RefGenomeDirectory '/genome.fasta'];
fastaFile = fastaread(fastaFilename);

nts='atcg';
rc='tagc';

[ChrStarts, glength, ~, ScafNames]= genomestats(RefGenomeDirectory);

%%
tic; fprintf(1,'Reading in gb file...\n ');

clear CDS
for i=1:length(ScafNames)
    
    ScafNames_i=ScafNames{i};
    
    fprintf(1,[ScafNames_i '...']);
    
    % GRAB HEADER FOR GB FILE
    f=find(ScafNames_i=='|',2,'last');
    if f > 1
        fn = ScafNames_i(f(1)+1:f(2)-1) ;
    else
        fn=ScafNames_i; % if header doesn't have NCBI formatted header
    end
    
    % FIND ANNOTATION SOURCE FOR SCAFFOLD SEQUENCE (GB or FASTA)
    gbfilename = [RefGenomeDirectory '/' fn '.gb'];
    
    if exist(gbfilename, 'file')
        
        % LOAD GENBANK
        genbankFile = genbankread(gbfilename);
        % check if sequence is empty
        if isfield(genbankFile,'Sequence')==0
            genbankFile.Sequence = lower(fastaFile(i).Sequence);
        end
        scafSeq = genbankFile.Sequence;
        
        
    else
        error(['Could not find a genebank file named' gbfilename]);
    end
   
    if isfield(genbankFile,'CDS')
        if length(genbankFile.CDS)>4;
            length(genbankFile.CDS)
            genes = locustag_from_text(genbankFile.CDS) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%     genes = div_add_in_nonCDS(genes, genbankFile.Features);
             CDS{i} = parse_all_locations_gb(genes, char([fastaFile.Sequence]+32)) ;  %also reverses strands in this, sorts by position
    %sort by position on genome
            [~,sortedpositions]=sort([CDS{i}.loc1]);
            CDS{i}=CDS{i}(sortedpositions);    
        else
            CDS{i} = {};
        end
    else
        CDS{i}={};
    end
    
    
    
end


%% save

save([RefGenomeDirectory '/cds_sorted'],'CDS')

return
