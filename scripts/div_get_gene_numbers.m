function genesN = div_get_gene_numbers(genestructure)

%Tami Lieberman
%May 2012
%Edited many times to deal with specific genomes
%Edited November 2016 to look for when the text is a number

genesN=zeros(length(genestructure),1);

tags={genestructure.locustag}';
tags
num_starts=cellfun(@locustag2num,tags); % <-- this is a key step see the function locustag2num
for i = 1:length(genestructure)
    if num_starts(i)>0
        %If than tag contains 'nc' that this is a noncoding element and
        %don't use it'
        if isempty(strfind(tags{i},'nc'));
            if tags{i}(end)>=48 & tags{i}(end)<58 %if the tag ends in a number
                genesN(i,1)=str2double(tags{i}(num_starts(i):end));
            else %the tag could end in a 'c' like it does in TB'
                genesN(i,1)=str2double(tags{i}(num_starts(i):end-1));
            end
        end
    end
end

%genesN=unique(genesN);  %not sure why this was here previously, commented
%out
%3/6/13

end