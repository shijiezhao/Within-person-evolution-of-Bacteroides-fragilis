function genelengths = get_gene_lengths(CDS)

allcds=CDS{:};
allgenenums=div_get_gene_numbers(CDS{1});

allcds=allcds(allgenenums>0); %usually RNA
allgenenums=allgenenums(allgenenums>0);


[allgenenums,codingsequencepostions]=sort(allgenenums);
allcds=allcds(codingsequencepostions);

for i=find(cellfun(@isempty,{allcds(:).loc2})>0)
    allcds(i).loc2=allcds(i).loc1+1000; %make arbitrary sizes for others
end

genelengths=zeros(1,max(allgenenums));
genelengths(allgenenums)=[allcds.loc2]-[allcds.loc1];

end