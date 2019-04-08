function IsolateTable = read_sample_table

%Line deliminter of .csv files must be \n\r


d = tdfreadunix('samples.csv',',') ;
%pwd
%d = csvread('samples.csv') ;


%disp(fieldnames(d))

if ~all(ismember({'Batch','Lane','Barcode','Sample','Alignments'},fieldnames(d)))
    error('Notice headers in samples.csv')
end

k = {} ;
for f=fieldnames(d)'
    if ischar(d.(f{1}))
        d.(f{1}) = cellstr(d.(f{1})) ;
    else
        d.(f{1}) = num2cell(d.(f{1})) ;
    end
    k(end+1,:) = {d.(f{1}){:}} ;
end
IsolateTable = cell2struct(k,fieldnames(d),1) ;

for i=1:length(IsolateTable)
    tmp = textscan(IsolateTable(i).Alignments,'%s ',1000) ;
    IsolateTable(i).Alignments = tmp{1} ;
    tmp = textscan(IsolateTable(i).Batch,'%s ',1000) ;
    IsolateTable(i).Batch = tmp{1} ;
end


end