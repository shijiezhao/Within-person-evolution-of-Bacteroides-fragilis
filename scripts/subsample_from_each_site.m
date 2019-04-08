function [SiteTimes, SiteFreqs] = subsample_from_each_site(isolateSites,observedMatrix,genotypeMatrix,...
    numSites, minsamples,  min_freq_for_assignment,max_remaining_room,max_to_ignore)

%This function does several things:
% 1) Removes locations with fewer than minsamples
% 2) Subsamples each location to the minimum of remaining samples
% 3) Assigns genotypes SiteTimes and SiteFreqs on subsampled data


%% subsample to make a smaller observed matrix, based on isolate_Sites

%decide how many to take from each Site
numPerSite=zeros(numSites,1);
for c=1:numSites
    numPerSite(c)=sum(isolateSites==c);
end
goodSites=find(numPerSite>=minsamples);
toTakeFromEachSite=min(numPerSite(goodSites));

%fprintf(1,[num2str(toTakeFromEachSite) ' ']);

%randomly choose samples to take
newIsolateSites=zeros(size(isolateSites));
for i=1:numel(goodSites)
    newIsolateSites(datasample(find(isolateSites==goodSites(i)),toTakeFromEachSite,'Replace',false))=goodSites(i);
end

%make new observed matrix
observedMatrix=observedMatrix(:,newIsolateSites>0);
newIsolateSites=newIsolateSites(newIsolateSites>0);

[SiteFreqs, SiteTimes] = assign_genotypes(observedMatrix,newIsolateSites,genotypeMatrix,numSites, min_freq_for_assignment,max_remaining_room,max_to_ignore, {});
