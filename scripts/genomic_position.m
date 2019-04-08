function gp = genomic_position(dfasta,pos)

gp=0.5*ones(size(pos)) ;


genestarts=[dfasta.loc1];
geneends=[dfasta.loc2];

for i=1:(numel(genestarts)-1)
    gp(genestarts(i):geneends(i))=i;
    gp(geneends(i)+1:genestarts(i+1)-1)=i+0.5;
end
gp(genestarts(i+1):geneends(i+1))=i+1;

return