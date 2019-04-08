function [siteFreqs, siteTimes, calls_for_phylogeny, names_for_phylogeny, genotypesByNumberMutations] = ...
    assign_genotypes(observedMatrix,isolateSites, genotypeMatrix,numsites,min_freq_for_assignment,max_remaining_room,max_to_ignore, SampleNames)


%% Assigns to genotypes to make siteFreqs and siteTimes


%Initialize
numpositions=size(genotypeMatrix,2);
mutsPerGenotype=sum(genotypeMatrix,2);
[~,genotypesByNumberMutations]=sort(mutsPerGenotype,'descend');

siteFreqs=zeros(length(genotypesByNumberMutations)+1,numsites); %last element is ancestral
siteTimes=zeros(length(genotypesByNumberMutations)+1,numsites); %last element is ancestral

%Optional stuff for plotting / troubleshooting
if ~isempty(SampleNames)
    calls_for_phylogeny=zeros(numpositions,1);
    names_for_phylogeny={};
    names_for_phylogeny{1}='Anc';
end

%For each each isolate
for i=1:numel(isolateSites)
    
    freqsleft=observedMatrix(:,i)';
    
    freqsfound=[];
    
    %try to assign to genotypes, starting with genotypes with most mutations
    for j=1:numel(genotypesByNumberMutations)
        
        genotypen=genotypesByNumberMutations(j);
        g=find(genotypeMatrix(genotypen,:)>0); %which mutations are in this genotype

        if min(freqsleft(g)) > min_freq_for_assignment
            
            genotypefreq=min(freqsleft(g));
            freqsleft(g)=freqsleft(g)-genotypefreq; %subtract this genotype from observation
            
            %record for data
            siteFreqs(j,isolateSites(i))=siteFreqs(j,isolateSites(i))+genotypefreq;
            siteTimes(j,isolateSites(i))=siteTimes(j,isolateSites(i))+1;
            freqsfound(end+1)=genotypefreq;
            
            if ~isempty(SampleNames)
                calls_for_phylogeny=[calls_for_phylogeny genotypeMatrix(genotypen,:)'];
                if sum(floor(freqsfound*100)==floor(genotypefreq*100))==0
                    names_for_phylogeny(end+1)={[SampleNames{i}(5:end) '-' num2str(floor(genotypefreq*100))]};
                else %another genotype already found in this sample with same frequency
                    names_for_phylogeny(end+1)={[SampleNames{i}(5:end) '-' num2str(floor(genotypefreq*100)) 'a' ]};
                end
            end
        end
    end
    
    %ancestral case
    if sum(freqsleft) < max_to_ignore & sum(freqsfound)< (1 - max_remaining_room)
        genotypefreq=1-sum(freqsfound);
        freqsfound(end+1)=genotypefreq;
        siteFreqs(end,isolateSites(i))=siteFreqs(end,isolateSites(i))+genotypefreq;
        siteTimes(end,isolateSites(i))=siteTimes(end,isolateSites(i))+1;
        
        if ~isempty(SampleNames)
            calls_for_phylogeny=[calls_for_phylogeny zeros(numpositions,1)];
            if sum(floor(freqsfound*100)==floor(genotypefreq*100))==0
                names_for_phylogeny(end+1)={[SampleNames{i}(5:end) '-' num2str(floor(genotypefreq*100))]};
            else %another genotype already found in this sample with same frequency
                names_for_phylogeny(end+1)={[SampleNames{i}(5:end) '-' num2str(floor(genotypefreq*100)) 'a' ]};
            end
        end
        
        
    end
    
end
