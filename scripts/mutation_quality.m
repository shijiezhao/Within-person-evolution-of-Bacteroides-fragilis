function [MutQual, MutQualIsolates] = mutation_quality(Call,Qual)



% MutQual -- Reports the minimum quality score in the best disagreement min(max for each type of call) 
% MutQualIsolates -- Reports the pair of isolates used to make this call


[Nmuts, NStrain] = size(Call) ;



MutQual = zeros(Nmuts,1) ; 
MutQualIsolates = zeros(Nmuts,2); 

for k=1:Nmuts
    
    if (length(unique([Call(k,:), 'N']))<=2) ; %if there is only one type of non-N call, skip this location
        MutQualIsolates(k,:) = 0; 
    else
        c1=Call(k,:) ; c2=c1' ;
        q2=Qual(k,:) ; q2=q1' ;
        g=c1~=c2 & c1~='N' & c2~='N' ; %find places where calls disagree
        positive_pos = find(g); 
        
        % get MutQual + logical index for where this occurred
        [MutQual(k), MutQualIndex] = max(min(q1(g),q2(g))) ;%min(q1(g),g2(g)) gives lower qual for each disagreeing pair of calls, we then find the best of these
        % store which strains were used to call this mutation 
        [strain_i, strain_j] = ind2sub(size(g), positive_pos(MutQualIndex));
        MutQualIsolates(k,:) = [strain_i, strain_j]; 
        
    end
end

MutQual(isnan(MutQual))=0;

return