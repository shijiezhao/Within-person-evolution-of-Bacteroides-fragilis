%% Step 1: Separate the samples by sampling time
TPS = {}; % Each entry in TPS is a collection of samples from a distinct time point
stool_ids=[]; % the corresponding stool sample ID
mk = 1;
for i = 1:length(SampleNames);
    tmp = strsplit(SampleNames{i},'_');
    if strcmp(tmp(1),Donor); % not from out_group
        st_id = str2num(tmp{2});
        if find(stool_ids==st_id);
            indx = find(stool_ids==st_id);
            TPS{indx} = [TPS{indx};i];
        else;
            stool_ids = [stool_ids;st_id];
            indx = find(stool_ids==st_id);
            TPS{indx} = [i];
        end
    end
end


%%
AFCt=[]; % Average fold-coverage for each time point
AFCt=[];
for ltp = 1:length(TPS);
    TempGoodsamples=TPS{ltp}';
    AFCt = [AFCt; mean(mean(all_coverage_per_bp(TempGoodsamples,:)))];
end
save('coverage.mat','AFCt')