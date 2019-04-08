%% 10. Make a tree for each SNP location
samplestoplot=samplestoplot(~in_outgroup);

if ~exist('tree_counting','dir')
    mkdir('tree_counting')
end

cd('tree_counting')

fid=fopen('for_tree_labeling.csv','w');
fprintf(fid,'chr,pos');
for i=1:numel(TreeSampleNames)
    fprintf(fid,[',' TreeSampleNames{i}]);
end
for i=1:numel(goodpos)
    fprintf(fid,['\n' num2str(contig_positions(positions_to_show_in_table(i),1)) ',' num2str(contig_positions(positions_to_show_in_table(i),2))]);
    for j=1:size(calls_for_tree,2)
        fprintf(fid,[',' calls_for_tree(i,j)]);
    end
end
fprintf(fid,'\n');
fclose(fid);
eval(['! python '  '../../../../scripts/countMutations.py ../' treename ' for_tree_labeling.csv'])
cd('..')