function IsolateTable = read_sample_names

d = tdfreadunix('sample_names.csv',',') ;

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

end