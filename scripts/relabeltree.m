function relabeltree(treefile,newtreefile,oldnames,newnames)


%%  Purpose: To relabel the tips of trees in a Newick-formatted file.
%Especially for using with phylip programs, which require tip labels which
%are <=10 characters long

%Tami Lieberman 2016

%%

fin=fopen(treefile,'r');
fout=fopen(newtreefile,'w');

x=fread(fin)';
newtree=x;

for i=1:numel(oldnames)
    j=strfind(char(newtree),[oldnames{i} ':']);
    oldnamelength=numel(oldnames{i});
    if ~isempty(j)
        if numel(j)==1 %& j < numel(x)-30
            newtree=[newtree(1:j-1) newnames{i} newtree((j+oldnamelength):end)];
        else
            error(['two instances of oldname: ' oldnames{i}]);
        end
    else
        fprintf(1,['Warning: No sample found in treefile with name :' oldnames{i} '\n']);
    end
end

fprintf(fout,newtree);

fclose(fin);
fclose(fout);
        
            
        

