function [mut_genenum, CDS, mut_annotations] = annotate_mutations_auto_gb(Positions,RefGenomeDirectory)


%% Written by Idan Yelin, Roy Kishony, Tami Lieberman, and Hattie Chung

% May need revisions for your reference genome



%% Intialize

% Load fasta 
fastaFilename = [RefGenomeDirectory '/genome.fasta'];
fastaFile = fastaread(fastaFilename); 

nts='atcg';
rc='tagc';

[~, ~, ~, ScafNames]= genomestats(RefGenomeDirectory);

mut_genenum = zeros(size(Positions,1),1);

%%
tic; fprintf(1,'annotate_mutations...\n ');



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
        if isempty(genbankFile.Sequence)
            genbankFile.Sequence = lower(fastaFile(i).Sequence); 
        end
        scafSeq = genbankFile.Sequence; 
        
        SOURCE_ANNOTATED = 1; 
        
    else
       error(['Could not find a genebank file named' gbfilename]);  
    end
    
    % IF GENBANK SOURCE, ANNOTATE
    
    if SOURCE_ANNOTATED
        genes = locustag_from_text(genbankFile.CDS) ;
        genes = div_add_in_nonCDS(genes, genbankFile.Features);
        CDS{i} = parse_all_locations_gb(genes, char([fastaFile.Sequence]+32)) ;  %also reverses strands in this, sorts by position
        
        %sort by position on genome
        [~,sortedpositions]=sort([CDS{i}.loc1]);
        CDS{i}=CDS{i}(sortedpositions);
        z = Positions(:,1)==i; % pulls positions that are on this SCAFFOLD
        mut_genenum(z) = genomic_position(CDS{i},Positions(z,2));
    else
        genesPlaceholder = make_annotation_placeholder(ScafNames_i, scafSeq); 
        CDS{i}=genesPlaceholder; 
        z = Positions(:,1)==i; 
        mut_genenum(z) = 1; 
    end
    
    
end



% Fill annotation for mutated positions 

mut_annotations=[] ;
for i=1:size(Positions,1)
    
    
    if ~mod(i,100), fprintf('\nOn %i\n', i); end
    
    mut_annotations(i).gene_num = mut_genenum(i) ;
    mut_annotations(i).scaffold = Positions(i,1) ;
    mut_annotations(i).pos = Positions(i,2) ;
    mut_annotations(i).ref=char(ref_sequences{mut_annotations(i).scaffold}(mut_annotations(i).pos)-32);
    
    nScf = Positions(i,1) ;
    if mut_genenum(i)==round(mut_genenum(i)) % intragenic
        cdf = CDS{nScf}(mut_genenum(i)) ;
        
        %if product is more than one line, reshape
        if size(cdf.product,1)>1
            cdf.product=reshape(cdf.product',size(cdf.product,1)*size(cdf.product,2),1)';
        end
        mut_annotations(i).gene       = cdf.gene ;
        mut_annotations(i).protein    = cdf.product ;
        mut_annotations(i).protein_id = cdf.protein_id ;
        mut_annotations(i).strand     = cdf.strand ; %1 indicates reverse strand, 0 forward strand
        mut_annotations(i).loc1       = cdf.loc1 ;
        mut_annotations(i).loc2       = cdf.loc2 ;
        mut_annotations(i).Sequence   = cdf.Sequence ;
        mut_annotations(i).note       = cdf.note ;
        mut_annotations(i).locustag   = cdf.locustag ;
        mut_annotations(i).translation = nt2aa(cdf.Sequence, 'GENETICCODE', 11,'ACGTOnly','F');
                
        
            
        if mut_annotations(i).strand
            p = double(mut_annotations(i).loc2) - double(Positions(i,2)) + 1;
        else
            p = double(Positions(i,2)) -double(mut_annotations(i).loc1) + 1;
        end
        mut_annotations(i).nt_pos = p ;
        
        aan = floor((p-1)/3) + 1 ;
        ncn = p-(aan-1)*3 ;
        codons=cell(4,1);
        AA='';
        mut_annotations(i).aa_pos=aan;
        
        if numel(mut_annotations(i).Sequence) >= aan*3 & mut_annotations(i).translation > 1;
            codon = mut_annotations(i).Sequence(aan*3-2:aan*3) ;
            for j=1:4 %for each nts
                if mut_annotations(i).strand
                    codon(ncn)=rc(j);
                else
                    codon(ncn)=nts(j);
                end
                codons{j}=codon;
                AA(j) = nt2aa(codon, 'ACGTOnly', false, 'ALTERNATIVESTARTCODONS','F', 'GENETICCODE', 11) ;
                %setting ACGTonly to false prevents nt2aa from throwing an error for non-ACTG calls
                %setting ALTERNATIVESTARTCODONS to false prevents nts from
                %being called as alternative start codons
            end
        end
        
        
        mut_annotations(i).codons   = codons ;
        mut_annotations(i).AA   = AA ;
        mut_annotations(i).NonSyn = length(unique(AA))>1 ;
    else %intergenic
        if floor(mut_genenum(i))>=1 % gene before position
            cdf = CDS{nScf}(floor(mut_genenum(i))) ;
            mut_annotations(i).gene1 = cdf.gene ;
            mut_annotations(i).locustag1 = cdf.locustag ;
            mut_annotations(i).distance1 = Positions(i,2) - cdf.loc2;
            mut_annotations(i).protein1    = cdf.product ;
            if cdf.strand==1
                mut_annotations(i).distance1 = mut_annotations(i).distance1 * -1;
            end
        end
        if ceil(mut_genenum(i))<=length(CDS{Positions(i,1)}) % gene after position
            cdf = CDS{nScf}(ceil(mut_genenum(i))) ;
            mut_annotations(i).gene2 = cdf.gene ;
            mut_annotations(i).locustag2 = cdf.locustag;
            mut_annotations(i).distance2 = cdf.loc1 - Positions(i,2);
            mut_annotations(i).protein2    = cdf.product ;
            if cdf.strand==0
                mut_annotations(i).distance2 = mut_annotations(i).distance2 * -1;
            end
        end
        mut_annotations(i).NonSyn = nan;
    end
    
end

return
