%% 8. Create a parsimony tree
% 8.1. Create matrix for plotting
samplestoplot=[1:numel(SampleNames(~in_outgroup))];
%samplestoplot=[1:numel(SampleNames)];
% samplestoplot=samplestoplot(~in_outgroup);
positions_for_tree = goodpos;
% positions_for_tree=sum(hasmutation,2)>0;
Calls_for_treei=majorNT(goodpos,samplestoplot);
% Calls_for_treei=Calls(goodpos,samplestoplot);
% outgroup_nts=ancnti(positions_for_tree);
% ogsamples=find(in_outgroup);
% outgroup_nts=mode(majorNT(positions_to_show_in_table,ogsamples),2);
Calls_for_treei = [ancnti(goodpos), Calls_for_treei];
calls_for_tree=zeros(size(Calls_for_treei));
calls_for_tree(Calls_for_treei>0)=NTs(Calls_for_treei(Calls_for_treei>0));
calls_for_tree(Calls_for_treei==0)='N';

if strcmp(Donor,'D44');
    Calls_for_treei(27,1)=1;
end
%outgroup_nts = extract_outgroup_mutation_positions(REFGENOMEFOLDER, p2chrpos(p(positions_for_tree),ChrStarts));
TreeSampleNames= {SampleNames{samplestoplot}};
% calls_for_tree = [outgroup_nts, calls_for_tree];
TreeSampleNames={'Reference' SampleNames{samplestoplot}};
for i=1:numel(TreeSampleNames)
    if length(TreeSampleNames{i}) > 10
        x=strsplit(TreeSampleNames{i},'_');
        TreeSampleNames(i)={[x{1}  x{2}]};
        disp(TreeSampleNames{i})
    end
    if TreeSampleNames{i}(1)=='B';
       TreeSampleNames(i)={['B2_' TreeSampleNames{i}(2:8)]}
    end
end

% 8.2 Plot parsimonious tree
mkdir('PhylogeneticTrees');
treename = generate_parsimony_tree(calls_for_tree, TreeSampleNames, TreefileName);
fprintf(1,'\nNow we are done with the tree\n');