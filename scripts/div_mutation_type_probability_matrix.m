function div_mutation_type_probability_matrix(GenomeDir,promotersize)


%Tami Lieberman May 2012, Updated May 2016
%Made from annotate_mutations_auto_gb and div_find_amino_acids

%Permutes every possible mutation on genome and returns matrices with the
%incidience of each type of mutations


genomewide_possibilities = zeros(4,4,5); %original nt, new nt, type
cds_possibilities = zeros(4,4,2); %original nt on coding strand, new nt on coding stand, type (N vs S only)

%types are, in order:  Nonsynonymous, Synonymous, Intergenic, Promoter, In poorly annotated gene

%%


[ChrStarts, glength, ~, ~, sequences]= genomestats(GenomeDir);

if ~exist([GenomeDir '/cds_sorted.mat'], 'file')
    read_gb_modified(REFGENOMEFOLDER)
end

load([GenomeDir '/cds_sorted'])



Positions=p2chrpos([1:glength]',ChrStarts);

nts='atcg';
rc='tagc';

rv=[2 1 4 3];

Sequences={};
for i=1:numel(sequences)
    if double(sequences{i}(1))<91
        [~,Sequences{i}]=ismember(double(sequences{i}),double('ATCG'));
    else
        [~,Sequences{i}]=ismember(double(sequences{i}),double(nts));
    end
end


%%

for i=1:size(Positions,1)

    ref=Sequences{Positions(i,1)}(Positions(i,2));
    
    
    if ref > 0
        
        nonref=1:4; nonref(ref)=[];
        
        if genenum(i) > 0 && genenum(i)==round(genenum(i)) % intragenic
            
            
            cdf = CDS{Positions(i,1)}(genenum(i)) ;
            
            if cdf.strand
                p = double(cdf.loc2) - double(Positions(i,2)) + 1;
                
            else
                p = double(Positions(i,2)) -double(cdf.loc1) + 1;
                if nts(ref) ~= cdf.Sequence(p)
                    disp(Positions(i))
                    disp(nts(ref))
                    disp(cdf.Sequence(p))
                end
            end
            
            
            aan = floor((p-1)/3) + 1 ;
            
            if numel(cdf.Sequence) >= aan*3
                
                ncn = p-(aan-1)*3 ;
                codon = cdf.Sequence(aan*3-2:aan*3) ;
                refAA= double(nt2aa(codon, 'ACGTOnly', false, 'ALTERNATIVESTARTCODONS','F', 'GENETICCODE', 11));
                for j=nonref %for each nts
                    
                    
                    %find amino acid
                    if cdf.strand
                        codon(ncn)=rc(j);
                        coi=rv(ref); %coding strand original index
                        cni=rv(j);%coding strand new index
                    else
                        codon(ncn)=nts(j);
                        coi=ref;
                        cni=j;
                    end
                    
                    %is same?
                    if double(nt2aa(codon, 'ACGTOnly', false, 'ALTERNATIVESTARTCODONS','F', 'GENETICCODE', 11)) ~= refAA
                        genomewide_possibilities(ref,j,1) = genomewide_possibilities(ref,j,1)+1;
                        cds_possibilities(coi,cni,1) = cds_possibilities(coi,cni,1)+1;
                    else
                        genomewide_possibilities(ref,j,2) = genomewide_possibilities(ref,j,2)+1;
                        cds_possibilities(coi,cni,2) = cds_possibilities(coi,cni,2)+1;
                        
                    end
                    %setting ACGTonly to false prevents nt2aa from throwing an error for non-ACTG calls
                    %setting ALTERNATIVESTARTCODONS to false prevents nts from
                    %being called as alternative start codons
                end
                
            else
                genomewide_possibilities(ref,nonref,5)=genomewide_possibilities(ref,nonref,5)+1;
            end
        else
            
            p=0;
            if floor(genenum(i))>=1 % gene before position
                cdf = CDS{Positions(i,1)}(floor(genenum(i))) ;
                distance1 = Positions(i,2) - cdf.loc2;
                if cdf.strand==1
                    distance1 =distance1 * -1;
                end
                if distance1 <  0 & distance1 > -1*promotersize
                    p=1;
                end
            end
            if ceil(genenum(i))<=length(CDS{Positions(i,1)}) % gene after position
                cdf = CDS{Positions(i,1)}(ceil(genenum(i))) ;
                distance2 = cdf.loc1 - Positions(i,2);
                if cdf.strand==0
                    distance2 = distance2 * -1;
                end
                if distance2 <  0 & distance2 > -1*promotersize
                    p=1;
                end
            end
            
            
            if p==1
                genomewide_possibilities(ref,nonref,4)=genomewide_possibilities(ref,nonref,4)+1;
            else
                genomewide_possibilities(ref,nonref,3)=genomewide_possibilities(ref,nonref,3)+1;
            end
        end
    end
end

probN = div_matrix2_6types(genomewide_possibilities(:,:,1)./ (genomewide_possibilities(:,:,1)+genomewide_possibilities(:,:,2)))/2;

%% 

save([GenomeDir '/dNdStools'], 'genomewide_possibilities', 'cds_possibilities', 'probN')