function ref_ID  = extract_outgroup_mutation_positions(ref_folder, positions)
    
    fastafile = [ref_folder '/genome.fasta']; 
    
    fr = fastaread(fastafile) ;
    
    ref_ID=zeros(length(positions),1);
    for i=1:length(positions)
        ref_ID(i) = fr(positions(i,1)).Sequence(positions(i,2)); 
    end
    
end