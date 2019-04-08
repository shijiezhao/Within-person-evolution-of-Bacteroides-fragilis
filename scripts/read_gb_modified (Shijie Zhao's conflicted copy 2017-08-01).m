function read_gb_modified(RefGenomeDirectory)


%% Written by Idan Yelin, Roy Kishony, Tami Lieberman, and Hattie Chung

% May need revisions for your reference genome



%% Intialize

% Load fasta
fastaFilename = [RefGenomeDirectory '/genome.fasta'];
fastaFile = fastaread(fastaFilename);

nts='atcg';
rc='tagc';

[ChrStarts, glength, ~, ScafNames]= genomestats(RefGenomeDirectory);

%%
tic; fprintf(1,'Reading in gb file...\n ');

clear CDS
for i=48:length(ScafNames)
    
    ScafNames_i=ScafNames{i};
    
    fprintf(1,[ScafNames_i '...\n']);
    
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
        if isempty(genbankFile.Sequence)
            genbankFile.Sequence = lower(fastaFile(i).Sequence);
        end
        scafSeq = genbankFile.Sequence;
        
        
    else
        error(['Could not find a genebank file named' gbfilename]);
    end
    
    genes = locustag_from_text(genbankFile.CDS) ;
    genes = div_add_in_nonCDS(genes, genbankFile.Features);
    CDS{i} = parse_all_locations_gb(genes, char([fastaFile.Sequence]+32)) ;  %also reverses strands in this, sorts by position
    
    %sort by position on genome
    if ~isempty(CDS{i})
        [~,sortedpositions]=sort([CDS{i}.loc1]);
        CDS{i}=CDS{i}(sortedpositions);
    end
    
    
end

Positions=p2chrpos([1:glength]',ChrStarts);

genenum = zeros(glength,1) ;

for i=1:numel(CDS)
    z = Positions(:,1)==i ;
    if ~isempty(CDS{i})
        genenum(z) = genomic_position(CDS{i},Positions(z,2)) ;
    else
        genenum(z)=0.5;
    end
end



%% save

save([RefGenomeDirectory '/cds_sorted'],'CDS', 'genenum')

return
