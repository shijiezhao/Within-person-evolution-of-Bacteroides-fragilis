function [genotypematrix, defined_genotype] = find_genotypes(fmatrix, avgcov, ismixedinfection,...
    strainmuts, min_snp_freq_single_sample, purity_cutoff_for_sorting, purity_cutoff_for_including_additional_singletons,...
    mincov,loose_majority_threshold, strict_majority_threshold, noise_buffer_for_singletons)



%initialization
defined_genotype=[]; %stores which sample defined each genotype, useful for troubleshooting
numgenotypes=0;
numpositions=size(fmatrix,1);
numdenovopositions=size(fmatrix,1);
genotypematrix=zeros(numpositions,numpositions); %each row is a genotype of a new strain, each column is a position
if ismixedinfection > 0
    % add in ancestor of each lineage
    numdenovopositions=numdenovopositions-(2*strainmuts);
    genotypematrix(1,:)=[zeros(1, numdenovopositions) zeros(1,strainmuts) ones(1,strainmuts)];
    genotypematrix(2,:)=[zeros(1, numdenovopositions) ones(1,strainmuts) zeros(1,strainmuts)];
    numgenotypes=2;
    defined_genotype(1:2)=0;
end


% sorting
forsorting=avgcov;
for i=1:numel(avgcov)
    %do pure strains first, followed by less pure strains
    %within each category (pure or less pure), sort in order of coverage
    freqs=fmatrix(1:numdenovopositions,i);
    if mean(freqs(freqs>(min_snp_freq_single_sample/2)))>purity_cutoff_for_sorting
        %min_snp_freq_single_sample is divided by two to make sure things close to cutoff don't result in pure call
        forsorting(i)=forsorting(i)+(max(avgcov)*2); %add a bonus for purity larger than max coverage
    end
end
[~,orderbycoverage]=sort(forsorting,'descend');





% going through each sample individually
for j=1:numel(orderbycoverage)
    
    i=orderbycoverage(j);
    
    if avgcov(i) > mincov
        
        genotype_strictly_defined=zeros(1,numpositions); genotype_loosely_defined=zeros(1,numpositions);
        genotype_strictly_defined(fmatrix(:,i)>strict_majority_threshold)=1;
        genotype_loosely_defined(fmatrix(:,i)>loose_majority_threshold)=1;
        
        
        %loose threshold -- min rather than mean because two related strains (one with additional mutation) might be in same sample
        majorstrainfreq=min(fmatrix(fmatrix(1:numdenovopositions,i)>=.5,i));
        
        if sum(genotype_strictly_defined)~=0  &  ...
                ~(sum(ismember(genotypematrix, genotype_strictly_defined, 'rows'))>0) & ...
                ~(sum(sum(abs(genotypematrix-repmat(genotype_loosely_defined,numpositions,1)),2)/sum(genotype_loosely_defined)<(1-majorstrainfreq))==1) & ...
                mean(fmatrix(genotype_loosely_defined>0,i)) > strict_majority_threshold
            %if at least one mutation in clear majority &
            %condition not found exactly before &
            %condition not found similar before, few of the mutations are on the edge of being in the the majority
            %really in strong majority
            %MAKE A NEW GENOTYPE
            numgenotypes=numgenotypes+1;
            defined_genotype(end+1)=i;
            genotypematrix(numgenotypes,:)=genotype_strictly_defined;
        elseif sum(fmatrix(1:numdenovopositions,i)>(min_snp_freq_single_sample-noise_buffer_for_singletons))==1 & ...
                sum(fmatrix(1:numdenovopositions,i)>min_snp_freq_single_sample)==1 & ...
                ~(sum(ismember(genotypematrix, fmatrix(:,i)'>min_snp_freq_single_sample, 'rows'))>0)
            %if mutation is really found alone and
            %this genotype not found before
            numgenotypes=numgenotypes+1;
            defined_genotype(end+1)=i;
            genotypematrix(numgenotypes,:)=fmatrix(:,i)>min_snp_freq_single_sample;
        end
        
        %if very strong majority strain and single other mutation at low frequency, add this minority genotype
        if (majorstrainfreq > purity_cutoff_for_including_additional_singletons & sum(fmatrix(:,i)>(1-majorstrainfreq+noise_buffer_for_singletons) & fmatrix(:,i)<majorstrainfreq)==1)
            genotype_strictly_defined(fmatrix(:,i) > (1-majorstrainfreq+noise_buffer_for_singletons))=1;
            if ~(sum(ismember(genotypematrix, genotype_strictly_defined, 'rows'))>0)
                %new genotype not found before
                numgenotypes=numgenotypes+1;
                defined_genotype(end+1)=i;
                genotypematrix(numgenotypes,:)=genotype_strictly_defined;
            end
        end
    end
    
end

genotypematrix=genotypematrix(1:numgenotypes,:);