function mutated = find_mutated_positions(Call)

%November 2012
%added variable strains so that an outgroup can be removed



mutated = zeros(size(Call,1),1) ; 

for k=1:size(Call,1);
    if (length(unique([Call(k,:), 0]))>2) ; %is there more than one type of non-N call after all the filters?
        mutated(k) = 1 ;
    end
end


return