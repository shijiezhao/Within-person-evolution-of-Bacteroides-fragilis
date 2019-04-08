function generate_phylip_input(calls, names, fileout)

if nargin < 3
    fid=1;
else
    fid = fopen(fileout,'w');
end

NStrain=size(calls,2);

%phylip header
fprintf(fid, ['   ' num2str(NStrain) '   ' num2str(size(calls,1)) '\n']);

%phylipinput=char(10+size(printCalls,1),size(printCalls,2));

for i=1:NStrain
    name=names{i}; 
    
    if numel(name) >= 10
        name=name(1:10);
    else
        for j=1:(10-numel(name))
            name=[name ' '];
        end
    end
    fprintf(fid, [name '   ' calls(:,i)' '\n']);
end