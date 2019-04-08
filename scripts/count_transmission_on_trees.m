function numtranmissions_from_lung = count_transmission_on_trees(genotypes, compartmenttimes, compartment_number_of_destination, ismultiplestraininfection, mintimes)

%Each row in genotypes is a genotype, each column is a mutation
%Each row in compartmenttimes is a genotype, each column is a compartment,
%values are the number of times it was found in that compartment

%Tami Lieberman, July 2015



%% variables you might want to change

source_compartment=1;


%% intialize
numtranmissions_from_lung=0;

genotypes_in_source=compartmenttimes(:,source_compartment)>=mintimes;
mutations_in_source=find(sum(genotypes(genotypes_in_source,:))>0);
    
genotypes_in_destination=compartmenttimes(:,compartment_number_of_destination)>=mintimes;
mutations_in_organ=find(sum(genotypes(genotypes_in_destination,:),1)>0);

unaccounted_for_shared_mutations=mutations_in_organ(ismember(mutations_in_organ,mutations_in_source));

genotypes_containing_shared_mutation=sum(genotypes(:,unaccounted_for_shared_mutations),2)>0 & (genotypes_in_destination | genotypes_in_source);


%% deal with special cases 

if ~ismultiplestraininfection
    %(1) ancestor of patient found in organ
    if sum(sum(genotypes(genotypes_in_destination,:),2)==0)>0 %ancestor in organ
        numtranmissions_from_lung=numtranmissions_from_lung+1;
   %(2) where the ancestor isn't found in organ but a unique descendent of it is
    elseif  sum(genotypes_in_destination & ~genotypes_containing_shared_mutation) > 0
        numtranmissions_from_lung=numtranmissions_from_lung+1;
    %(3) no genotype, including ancestor found more than
    %MIN_OBSERVATIONS_FOR_TRANSMISSION_ANALYSIS -- infer a single
    %transmission event of ancestor
    % (alternatively this could have been chosen to be an excluded case)
    elseif sum(genotypes_in_destination) == 0 
        numtranmissions_from_lung=numtranmissions_from_lung+1;
    end
end


%% iteratively account for mutations

while numel(unaccounted_for_shared_mutations)>0
    numtranmissions_from_lung=numtranmissions_from_lung+1;
    %find genotypes that contribute to these
    %unshared mutations
    %pick genotypes with fewest mutations first
    %first look for genotypes in organ 
    
    genotypes_containing_shared_mutation=find(sum(genotypes(:,unaccounted_for_shared_mutations),2)>0 & (genotypes_in_destination));
        
    if isempty(genotypes_containing_shared_mutation)
        error(1,'did not find a genotype with unaccoutned for mutation in destination organ');
    end
    
    minmuts=min(sum(genotypes(genotypes_containing_shared_mutation,:),2));
    mingenotype=genotypes_containing_shared_mutation(find(sum(genotypes(genotypes_containing_shared_mutation,:),2)==minmuts,1));
 
    %remove all unaccounted for mutations carried by this minimal genotype
    unaccounted_for_shared_mutations(ismember(unaccounted_for_shared_mutations,find(genotypes(mingenotype,:))))=[];
    
    
end




%% sanity check

if numtranmissions_from_lung==0
    error('no transmissions')
end
