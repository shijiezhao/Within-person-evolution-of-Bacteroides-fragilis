function maxdist = find_max_distance_between_genotypes(genotypeMatrix)

if isempty(genotypeMatrix)
    maxdist=0;
else
    maxdist=1;
    if size(genotypeMatrix,1) > 1
        for i=1:size(genotypeMatrix,1)-1
            for j=(i+1):size(genotypeMatrix,1)
                if sum(genotypeMatrix(i,:)~=genotypeMatrix(j,:)) > maxdist
                    maxdist = sum(genotypeMatrix(i,:)~=genotypeMatrix(j,:));
                end
            end
        end
    end
end

end