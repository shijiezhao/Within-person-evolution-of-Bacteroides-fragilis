function gp = find_genomic_position(dfasta,pos)

gp=zeros(size(pos)) ;

for i=1:length(pos)
    
    f1 = find(pos(i)>=[dfasta.loc1],1,'last') ;
    
    if isempty(f1)
        gp(i) = 0.5 ;
    elseif pos(i)<=dfasta(f1).loc2
        gp(i) = f1 ;
    else
        gp(i) = f1+0.5 ;
    end
    
end

return