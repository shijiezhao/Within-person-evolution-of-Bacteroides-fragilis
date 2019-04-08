function [ltagnumbers, intergenic] =get_locustags(p,codingsequences)

cgenes=genomic_position(codingsequences{1} ,p);

%remove intragenic regions
intergenic=cgenes~=floor(cgenes);
intragenic=find(~intergenic);
cgenes(intergenic)=[];

ltagnumbers=div_get_gene_numbers(codingsequences{1}(cgenes));

intergenic(intragenic(ltagnumbers==0))=1;  %mostly for RNA
ltagnumbers=ltagnumbers(ltagnumbers>0);


end

