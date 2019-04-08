function clickable_snp_table_no_annotations(contig_positions, callsi, cnts, SampleNames, ScafNames, MutQual, QualSort)

NTs='ATCG';

calls=repmat('N',size(callsi));
calls(callsi>0)=char(NTs(callsi(callsi>0)));

nump=size(cnts,2);

%generate table


%set column names
colnames={'Qual','Contig','Pos'};
nonsamplecols=numel(colnames);
for i=1:numel(SampleNames)
    colnames{end+1}=SampleNames{i};
end

%set column widths
widths=num2cell([48, 48, 100, ones(1,numel(SampleNames))*26]);



tabledata=cell(nump,numel(colnames));


for i=1:nump
    % ___ generate table ___ %
    
    tabledata(i,1:nonsamplecols)=[{MutQual(i)} contig_positions(i,1) contig_positions(i,2)];
    for j = 1:numel(SampleNames)
        tabledata(i,nonsamplecols+j) = cellstr(calls(i,j));
    end
end


positions = 1:size(tabledata,1);
oldtable=tabledata;
if QualSort==1
    [tabledata, positions] = sortrows(oldtable,1);
end

%display table
figure();clf;hold on;
set(gcf, 'Position',[10         50        1250         550]);
uicontrol('Style','text','Position',[400 45 120 20],'String','Vertical Exaggeration')
t = uitable('Units','normalized','Position',[0 0 1 .97], 'Data', tabledata,...
    'ColumnName', colnames,...
    'RowName',[], ...
    'CellSelectionCallback',@mut_matix_clicked, ...
    'ColumnWidth',widths);

fid=fopen('snptable.csv','w');
fprintf(fid,'Contig,Pos,');
for i=1:numel(SampleNames)
    fprintf(fid,[SampleNames{i} ',']);
end
fprintf(fid,'\n');
for i=1:size(contig_positions,1)
    fprintf(fid,[num2str(contig_positions(i,1)) ',' num2str(contig_positions(i,2)) ',']);
    for j=1:numel(SampleNames)
        if callsi(i,j)>0
            fprintf(fid,[NTs(callsi(i,j)) ',']);
        else
            fprintf(fid,'N,');
        end
    end
    fprintf(fid,'\n');
end
fclose(fid);


%intialize divbarchartswindow
figure(660);clf;

    function mut_matix_clicked(~, event)
                
        rc = event.Indices ;
        
        ind=positions(rc(1));

        %Bar charts of counts

        
        div_bar_charts(squeeze(cnts(:,ind,:)), SampleNames)
        disp(ind)

    end


end